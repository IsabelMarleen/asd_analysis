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
gene_info <- read.delim( file.path( path, "genes.tsv" ), header=FALSE, as.is=TRUE ) %>%
  mutate(unique = case_when(
  duplicated(V2) | duplicated(V2, fromLast=T) ~ paste(V2, V1, sep="_"),
  TRUE ~ V2))
rownames(counts) <- gene_info$unique
colnames(counts) <- readLines( file.path( path, "barcodes.tsv" ) )


sampleTable <-
  cellinfo %>% select( sample : RNA.Integrity.Number ) %>% unique
sampleTable

# extracting gene expression is much faster in column-sparse format:
Tcounts <- as(t(counts), "dgCMatrix")   #  fast:   Tcounts[, "SYN1"]

Ccounts <- as(counts, "dgCMatrix")      #  fast:   Ccounts[, 1337]   and  colSums(Ccounts)



# find informative genes    (rsession goes up to 25 GB RAM [htop])
sfs <- colSums(counts)
norm_counts <- t(t(Ccounts) / colSums(Ccounts))
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
umap_euc <- uwot::umap( pca$x, spread = 10, n_threads = 40)

umap_cos <- uwot::umap( pca$x, metric = "cosine", spread = 10, n_threads = 40)






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
nn_cells <- FNN::get.knn( pca$x, k = 50)
adj <- Matrix(0, nrow = nrow(pca$x), ncol = nrow(pca$x)) # has to be sparse, otherwise takes 80 GB of RAM
for(i in 1:ncol(nn_cells$nn.index))
  adj[ cbind(1:nrow(pca$x), nn_cells$nn.index[, i]) ] <- 1
for(i in 1:ncol(nn_cells$nn.index))
  adj[ cbind(nn_cells$nn.index[, i], 1:nrow(pca$x)) ] <- 1
cl_louvain <- cluster_louvain(  graph_from_adjacency_matrix(adj, mode = "undirected") )

tmp_clusters <- cl_louvain$membership
tmp_clusters <- case_when(tmp_clusters %in% c(5, 3, 2, 8, 18, 20, 23, 19, 17) ~ 5, TRUE ~ tmp_clusters)




# plot clusters
df2 <- data.frame(cell = colnames(counts),
                  umap_euc,
                  cluster = factor(tmp_clusters),
                  stringsAsFactors = F) %>%
  left_join( select(cellinfo, cell, cluster) %>% rename(paper_cluster = cluster), by="cell" )
labels_louvain2 <- df2 %>% group_by(cluster) %>% summarise(X1 = mean(X1), X2=mean(X2))
labels_paper2 <- df2 %>% group_by(paper_cluster) %>% summarise(X1 = mean(X1), X2=mean(X2))
p_louv <- ggplot() + 
  geom_point(data=df2,
             aes(X1, X2, col = cluster),
             size = .05) + coord_fixed()+
  geom_label(data = labels_louvain2,
             aes(X1, X2, col = cluster, label = cluster),
             fontface = "bold")+ theme(legend.position = "none")
p_papcl <- ggplot() + 
  geom_point(data=df2,
             aes(X1, X2, col = paper_cluster),
             size = .05) + coord_fixed()+
  geom_label(data = labels_paper2,
             aes(X1, X2, col = paper_cluster, label = paper_cluster),
             fontface = "bold")+ theme(legend.position = "none")
p_louv # clusters from louvain
p_papcl # clusters from paper





# clean out doublets and ambiguous cells ----------------------------------


# number of NN from different cluster:
nn_inothercluster <- colSums(
  matrix(tmp_clusters[ t(nn_cells$nn.index) ],
         ncol = nrow(nn_cells$nn.index))   != 
  matrix(rep(tmp_clusters, each = ncol(nn_cells$nn.index)),
         ncol = nrow(nn_cells$nn.index)) )
 

# in silico doublets: randomly draw cells from different clusters and pool their UMIs to form a "synthetic" doublet:
cellsA <- c()
cellsB <- c()
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



# percentage of synthetic doublets doublets in neighborhood for each cell:
dblts_perc <- rowMeans( nn > nrow(pca$x) )[ 1:nrow(pca$x) ]



# Run UMAP  with Annoy's output
ump2 <- uwot::umap( NULL, nn_method = list( idx=nn, dist=nn_dists), 
                    n_threads=40, spread = 15, verbose=TRUE )



is_synth <- 1:nrow(ump2) > nrow(pca$x)






  



# Cluster contributions to doublets ---------------------------------------

gg <- ggplot_build(p_louv)
cl_cols <- unique(gg$data[[2]][c("label","colour")])

plot_grid(plotlist = c(list(p_louv),
lapply(c(21,1, 9, 6, 13), function(cl){
  ggplot()+coord_fixed()+
  geom_point(data=data.frame(ump2[1:nrow(pca$x),]), aes(X1, X2), col="grey", size=.05)+
  geom_point(data = data.frame(ump2[is_synth,],
                               cl_contributed = tmp_clusters[cellsA] == cl | tmp_clusters[cellsB] == cl),
             aes(X1, X2, col = cl_contributed), size=.05) +
   scale_color_manual(values = c(`FALSE`="grey", `TRUE`=cl_cols[cl_cols$label==cl,"colour"])) +
    
    geom_label(data = labels_df2[labels_df2$cluster == cl,],
               aes(X1, X2, label = cluster, fill=cluster),
               fontface = "bold")+ theme(legend.position = "none") +
   scale_fill_manual(values = cl_cols[cl_cols$label==cl,"colour"]) +
  
   ggtitle(paste0("Doublets with contribution from cluster ", cl)) + theme(legend.position = "none")
})
 ))



# synthetic doublets have synthetic doublets in their neighborhood
a <- FNN::get.knn(rbind(pca$x, doublet_pcs), k = 50)
data.frame(perc_dbl = rowMeans( a$nn.index > nrow(pca$x) ), is_synth = 1:nrow(a$nn.index) > nrow(pca$x)) %>%
  ggplot() + geom_histogram(aes(perc_dbl, fill = is_synth), alpha=.3) + coord_cartesian(ylim = c(0, 14000))+
  facet_wrap(~is_synth)


plot(rowMeans( a$nn.index > nrow(pca$x) ), pch=20, cex=.4); abline(h = nrow(pca$x))










# DESeq -------------------------------------------------------------------

# cells for which to compare ASD vs control:
sel <- tmp_clusters==5  & dblts_perc == 0  & nn_inothercluster <= 1



pseudobulks <- as.matrix(t( fac2sparse(cellinfo$sample[sel]) %*% t(Ccounts[, sel]) ))

coldat <- filter(sampleTable, sample %in% colnames(pseudobulks)) %>% 
  mutate(individual = factor(individual),
         diagnosis = factor(diagnosis, levels = c("Control", "ASD")),
         region    = factor(region))
rownames(coldat) <- coldat$sample


library(DESeq2)
library(BiocParallel)


dds <- DESeqDataSetFromMatrix( pseudobulks,
                               coldat[colnames(pseudobulks), ],
                               design = ~  diagnosis )
# I tested that we do not need interactions between sex, region and diagnosis with
# DESeq's LTR.
dds <- DESeq(dds, 
             parallel=TRUE, BPPARAM=MulticoreParam(20))
resultsNames(dds)
results(dds, name = "diagnosis_ASD_vs_Control") %>% as.data.frame() %>% rownames_to_column("Gene") %>%
  filter(padj < .1) %>% arrange(desc(log2FoldChange)) %>% dim


plot(
  -log10(results(dds_32364, name = "diagnosis_ASD_vs_Control")$padj),
  -log10(results(dds_23409, name = "diagnosis_ASD_vs_Control")$padj),
  pch=20, cex=.5, asp=1);abline(v=1, h=1, lty=2); abline(0,1)

# It looks like we are gaining a lot of power by removing ambiguous cells!





plotCounts(dds, "TTTY10", intgroup = c("sex", "region", "diagnosis"))
g <- "GADD45G"
data.frame(umap_euc, cellinfo, Gene = Tcounts[, g], sfs, sel) %>%
   filter(sel) %>%
  ggplot(aes(X1, X2, col = Gene / sfs / mean(1/sfs))) + geom_point(size=.1)+coord_fixed()+
  col_pwr_trans(1/2, g) + facet_wrap(~ region + diagnosis)





degs <- results(dds, name = "diagnosis_ASD_vs_Control") %>% as.data.frame() %>% rownames_to_column("Gene") %>%
  filter(padj < .1) %>% arrange(desc(log2FoldChange)) %>% pull(Gene)


deg_pca <- irlba::prcomp_irlba( x = sqrt(t(norm_counts[degs, sel])),
                            n = 20,
                            scale. = TRUE)
deg_umap  <- uwot::umap( deg_pca$x, n_neighbors=30, spread = 5, n_threads = 40, verbose = TRUE)

g <- "WNT7B"
data.frame(deg_umap, cellinfo[sel, ], Gene = Tcounts[sel, g], sfs=sfs[sel]) %>%
  ggplot(aes(X1, X2, col = diagnosis))+geom_point(size=1) + coord_fixed()
  ggplot(aes(X1, X2, col = Gene/sfs/mean(1/sfs)))+geom_point(size=1) + coord_fixed()+
    col_pwr_trans(1/2, g)









data.frame(umap_euc,
           sel
           ) %>% ggplot(aes(X1, X2, col = sel))+geom_point(size=.1)+coord_fixed()

model.matrix(~ sex + region + diagnosis, data = coldat[colnames(pseudobulks), ])


















