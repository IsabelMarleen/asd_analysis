# This script crashes with 16 GB or less of RAM.
# Rsession will use 30 GB or RAM in the long-run, not sure about peaks.

library( tidyverse )
library( Matrix )
library( irlba )
library( uwot )
library( FNN )
library( igraph )
library( cowplot )

# convenience functions such as col_pwr_trans, rowVars_spm etc.:
scr_dir <- "/home/frauhammer/sc_methods_dev/src/"
source(file.path(scr_dir, "functions_universal.R"))


path <- "/home/frauhammer/sds_copy/ASD/"
cellinfo <- read.delim( file.path( path, "meta.txt" ), stringsAsFactors=FALSE )
counts <- readMM( file.path( path, "matrix.mtx" ) )
rownames(counts) <- read.delim( file.path( path, "genes.tsv" ), header=FALSE )$V2
colnames(counts) <- readLines( file.path( path, "barcodes.tsv" ) )


sampleTable <-
  cellinfo %>% select( sample : RNA.Integrity.Number ) %>% unique
sampleTable

# extracting gene expression is much faster in column-sparse format:
Tcounts <- as(t(counts), "dgCMatrix")   #  fast:   Tcounts[, "SYN1"]

Ccounts <- as(counts, "dgCMatrix")      #  fast:   Ccounts[, 1337]






# find informative genes    (rsession goes up to 25 GB RAM [htop])
sfs <- colSums(counts)
norm_counts <- t(t(counts) / colSums(counts))
poisson_vmr <- mean(1/sfs)
gene_means <- rowMeans( norm_counts )
gene_vars <- rowVars_spm( norm_counts )
cells_expressing <- rowSums( counts != 0 )
is_informative <- gene_vars/gene_means > 1.5 * poisson_vmr  &  cells_expressing > 100

plot(gene_means, gene_vars/gene_means, pch=".", log = "xy")
points(gene_means[is_informative], (gene_vars/gene_means)[is_informative], pch=".", col = "red" )



# PCA and UMAP   (goes up to 18 GB of RAM)
pca <- irlba::prcomp_irlba( x = sqrt(t(norm_counts[is_informative,])),
                            n = 40,
                            scale. = TRUE)
umap_euc <- uwot::umap( pca$x, spread = 2, n_threads = 40)

umap_cos <- uwot::umap( pca$x, metric = "cosine", spread = 2, n_threads = 40)



# Clusters from paper -----------------------------------------------------

df <- data.frame(cell = colnames(norm_counts), umap_euc, stringsAsFactors = F) %>%
  left_join(cellinfo, by = "cell") %>%
  mutate(individual = factor(individual))
# some clusters from the paper are mixing in UMAP, i.e. are not robust to preprocessing, so we merge them:
df <- df %>% mutate(cluster = case_when(
  grepl("^Neu-NRGN-", cluster) ~ "Neu_NRGN",
  grepl("^AST-", cluster) ~ "AST",
  TRUE ~ cluster
))

labels_df <- df %>% group_by(cluster) %>% summarise(X1 = mean(X1), X2=mean(X2))

p_cl <-   ggplot() + 
  geom_point(data=df,
             aes(X1, X2, col = cluster),
             size = .05) + coord_fixed()+
  geom_label(data = labels_df,
             aes(X1, X2, col = cluster, label = cluster),
             fontface = "bold")+ theme(legend.position = "none")
p_cl





# Clusters and doublets ---------------------------------------------------


g <- "SYT1"    # Neuron marker
g <- "SLC1A2"  # Astrocyte marker

ggplot()+
  geom_point(data = data.frame(umap_euc, Gene = Tcounts[, g], sfs = sfs),
             aes(X1, X2, col = Gene / sfs / mean(1/sfs)), size=.1) +
  coord_fixed() + col_pwr_trans(1/2, g) +
  geom_label(data = labels_df,
             aes(X1, X2, label = cluster),
             fontface = "bold")




# Louvain clusters
nn <- FNN::get.knn( pca$x, k = 50)
adj <- Matrix(0, nrow = nrow(pca$x), ncol = nrow(pca$x)) # has to be sparse, otherwise takes 80 GB of RAM
for(i in 1:ncol(nn$nn.index))
  adj[ cbind(1:nrow(pca$x), nn$nn.index[, i]) ] <- 1
for(i in 1:ncol(nn$nn.index))
  adj[ cbind(nn$nn.index[, i], 1:nrow(pca$x)) ] <- 1
cl_louvain <- cluster_louvain(  graph_from_adjacency_matrix(adj, mode = "undirected") )







# in silico doublets: randomly draw cells from different clusters and pool their UMIs to form a "synthetic" doublet:
cellsA <- sample(1:ncol(counts), 5000)
cellsB <- sample(1:ncol(counts), 5000)

cellsA <- c()
cellsB <- c()
tmp_clusters <- cl_louvain$membership
tmp_clusters <- case_when(tmp_clusters %in% c(5, 3, 2, 8, 18, 20, 23, 19, 17) ~ 5, TRUE ~ tmp_clusters)
for(i in 1:10000){
  id1 <- sample(1:ncol(counts), 1)
  id2 <- sample( which(  
    # sample amongst cells from same sample and different cluster:
    cellinfo$sample == cellinfo$sample[id1]  &  tmp_clusters != tmp_clusters[id1]
    ),1)
  cellsA <<- c(cellsA, id1)
  cellsB <<- c(cellsB, id2)
}



doublet_raw <- Ccounts[, cellsA] + Ccounts[, cellsB]
doublet_pcs <- predict(pca,
                       newdata = sqrt( (t(doublet_raw) / colSums(doublet_raw))[, is_informative] ))



# Alternative 1 (clearer):
a <- FNN::get.knn(rbind(pca$x, doublet_pcs), k = 50)
nn <- a$nn.index
nn_dists <- a$nn.dist

# Alternative 2 (faster):
library( RcppAnnoy )
featureMatrix <- rbind(pca$x, doublet_pcs); k_nn <- 50
annoy <- new( AnnoyEuclidean, ncol(featureMatrix) )
for( i in 1:nrow(featureMatrix) ) 
  annoy$addItem( i-1, featureMatrix[i,] )
annoy$build( 50 ) # builds a forest  of n_trees trees. More trees gives higher precision when querying.
nn <- t( sapply( 1:annoy$getNItems(), function(i) annoy$getNNsByItem( i-1, k_nn) + 1 ) )
nn_dists <- sapply( 1:ncol(nn), function(j) sqrt( rowSums( ( featureMatrix - featureMatrix[ nn[,j], ] )^2 ) ) )


# Run UMAP  with Annoy's output
ump2 <- uwot::umap( NULL, nn_method = list( idx=nn, dist=nn_dists), 
                    n_threads=40, spread = 15, verbose=TRUE )



# plots
p_g1 <- data.frame(ump2[1:nrow(pca$x),], Gene = Tcounts[, "SYT1"]/sfs/mean(1/sfs)) %>% 
  ggplot(aes(X1, X2, col = Gene))+geom_point(size=.05)+coord_fixed()+col_pwr_trans(1/2, "Gene")

p_synth <- ggplot() + coord_fixed()+
  geom_point(data = data.frame(ump2[1:nrow(pca$x),]),
             aes(X1, X2),
             size = .05, col = "grey")+
  geom_point(data = data.frame(ump2[(nrow(pca$x)+1):nrow(ump2),]),
             aes(X1, X2),
             size = .05, col = "red") + ggtitle("10000 artificial doublets")



tmp_df <- data.frame(ump2[1:nrow(pca$x), ], sfs, cellinfo$cluster) 
tmp_ldf <- tmp_df %>% group_by(cellinfo.cluster) %>% summarise(X1 = mean(X1), X2=mean(X2))
p_umi <- ggplot()+coord_fixed()+geom_point(data=tmp_df, aes(X1, X2, col = log10(sfs)), size=.05)+
  scale_color_gradientn(name = "log10_nUMI", colours = rev(rje::cubeHelix(100))[5:100])+
  geom_label(data = tmp_ldf, aes(X1, X2, label = cellinfo.cluster)) + ggtitle("Neurons have high transcription")
p_umi

dblts_perc <- rowMeans( nn > nrow(pca$x) )[ 1:nrow(pca$x) ]
p_dbl <- ggplot() + coord_fixed()+
  geom_point(data = data.frame(ump2[1:nrow(pca$x),]),
             aes(X1, X2),
             size = .05, col = "grey")+
  geom_point(data = data.frame(ump2[1:nrow(pca$x),][dblts_perc > .05,]),
             aes(X1, X2),
             size = .05, col = "black")+
  ggtitle("Cells with synth. doublets as NN")
p_dbl



# plot clusters
df2 <- data.frame(ump2[1:nrow(pca$x),], cluster = factor(tmp_clusters))
labels_df2 <- df2 %>% group_by(cluster) %>% summarise(X1 = mean(X1), X2=mean(X2))
p_c <- ggplot() + 
  geom_point(data=df2,
             aes(X1, X2, col = cluster),
             size = .05) + coord_fixed()+
  geom_label(data = labels_df2,
             aes(X1, X2, col = cluster, label = cluster),
             fontface = "bold")+ theme(legend.position = "none")
p_c




plot_grid(p_c, p_synth, p_dbl, p_umi, ncol=2)









is_synth <- 1:nrow(ump2) > nrow(pca$x)

p_5 <- ggplot()+coord_fixed()+
  geom_point(data=data.frame(ump2[1:nrow(pca$x),]), aes(X1, X2), col="grey", size=.05)+
  geom_point(data = data.frame(ump2[is_synth,],
                               has_5 = tmp_clusters[cellsA] == 5 | tmp_clusters[cellsB] == 5),
             aes(X1, X2, col = has_5), size=.05)+ggtitle("Synthetic doublets with cell from cl5")

p_large <- ggplot()+coord_fixed()+
  geom_point(data=data.frame(ump2[1:nrow(pca$x),]), aes(X1, X2), col="grey", size=.05)+
  geom_point(data = data.frame(ump2[is_synth,],
                               has_large = sfs[cellsA] > 10^3.65 | sfs[cellsB] > 10^3.65),
             aes(X1, X2, col = has_large), size=.05) + ggtitle("Synthetic doublets containing large cells")

plot_grid(p_c+ggtitle("Louvain Clusters used to generate synth. doublets"),
          p_umi,
          p_large,
          ncol=3)




