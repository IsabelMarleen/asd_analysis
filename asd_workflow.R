# This script crashes with 16 GB or less of RAM.
# Rsession will use 30 GB or RAM in the long-run, not sure about peaks.





# Load --------------------------------------------------------------------


library( tidyverse )
library( Matrix )
library( irlba )
library( uwot )
library( FNN )
library( igraph )
library( cowplot )

# rowVars for sparse matrices:
colVars_spm <- function( spm ) {
  stopifnot( is( spm, "dgCMatrix" ) )
  ans <- sapply( seq.int(spm@Dim[2]), function(j) {
    mean <- sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
    sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) +
      mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) } ) / ( spm@Dim[1] - 1 )
  names(ans) <- spm@Dimnames[[2]]
  ans
}
rowVars_spm <- function( spm ) {
  colVars_spm( t(spm) )
}


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







# Preprocessing -----------------------------------------------------------

# load (or re-execute everything in this section):
sfs <- colSums(counts)
norm_counts <- t(t(Ccounts) / colSums(Ccounts))
rownames(norm_counts) <- rownames(Ccounts)
load(file.path("~", "asd_analysis", "savepoint", "umap_euc_spread10.RData"))




# informative genes, PCA, UMAP:
poisson_vmr <- mean(1/sfs)
gene_means <- rowMeans( norm_counts )
gene_vars <- rowVars_spm( norm_counts )
cells_expressing <- rowSums( counts != 0 )
is_informative <- gene_vars/gene_means > 1.5 * poisson_vmr  &  cells_expressing > 100
plot(gene_means, gene_vars/gene_means, pch=".", log = "xy")
points(gene_means[is_informative], (gene_vars/gene_means)[is_informative], pch=".", col = "red" )
pca <- irlba::prcomp_irlba( x = sqrt(t(norm_counts[is_informative,])),
                            n = 40,
                            scale. = TRUE)
umap_euc <- uwot::umap( pca$x, spread = 10, n_threads = 40)
umap_cos <- uwot::umap( pca$x, metric = "cosine", spread = 10, n_threads = 40)

# save(umap_euc,
#      file = file.path("~", "asd_analysis", "savepoint", "umap_euc_spread10.RData"))






# Clusters ---------------------------------------------------


# load (or re-execute everything in this section):
load(file.path("~", "asd_analysis", "savepoint", "clusters.RData"))






# find NN for each cell:
library( RcppAnnoy )
featureMatrix <- pca$x; k_nn <- 50
annoy <- new( AnnoyEuclidean, ncol(featureMatrix) )
for( i in 1:nrow(featureMatrix) ) 
  annoy$addItem( i-1, featureMatrix[i,] )
annoy$build( 50 ) # builds a forest  of n_trees trees. More trees gives higher precision when querying.
nn_cells <- t( sapply( 1:annoy$getNItems(), function(i) annoy$getNNsByItem( i-1, k_nn) + 1 ) )
nndists_cells <- sapply( 1:ncol(nn_cells), function(j) sqrt( rowSums( ( featureMatrix - featureMatrix[ nn_cells[,j], ] )^2 ) ) )
rm(featureMatrix, annoy)

# cluster on nearest neighbor graph (Louvain):
adj <- Matrix(0, nrow = nrow(pca$x), ncol = nrow(pca$x)) # has to be sparse, otherwise takes 80 GB of RAM
for(i in 1:ncol(nn_cells))
  adj[ cbind(1:nrow(pca$x), nn_cells[, i]) ] <- 1
for(i in 1:ncol(nn_cells))
  adj[ cbind(nn_cells[, i], 1:nrow(pca$x)) ] <- 1
cl_louvain <- cluster_louvain(  graph_from_adjacency_matrix(adj, mode = "undirected") )

# merge clusters that are separated by patient heterogeneity:
tmp_clusters <- cl_louvain$membership
tmp_clusters <- case_when(tmp_clusters %in% c(5, 6, 8, 1, 10, 20, 2) ~ 5, TRUE ~ tmp_clusters) # excit. Ns
tmp_clusters <- case_when(tmp_clusters %in% c(11, 15, 19) ~ 11, TRUE ~ tmp_clusters) # astrocytes
tmp_clusters <- case_when(tmp_clusters %in% c(3, 9, 18) ~ 3, TRUE ~ tmp_clusters) # OPCs





# Louvain clusters 
p_louv <- ggplot()+ coord_fixed() +
  geom_point(data = data.frame(umap_euc, cl=factor(tmp_clusters)),
             aes(X1, X2, col = cl), size = .1) +
  geom_label(data = group_by(data.frame(umap_euc, cl=factor(tmp_clusters)), cl) %>%summarise(X1=mean(X1), X2=mean(X2)), 
             aes(X1, X2, label = cl))
p_louv

# clusters from paper
p_paper <- ggplot()+ coord_fixed()+
  geom_point(data =data.frame(cell = colnames(counts), umap_euc) %>%
               left_join(select(cellinfo, cell, cluster), by="cell"),
             aes(X1, X2, col = cluster), size = .1) +
  geom_label(data = data.frame(cell = colnames(counts), umap_euc) %>%
               left_join(select(cellinfo, cell, cluster), by = "cell") %>% group_by(cluster) %>%
               summarise(X1=mean(X1), X2=mean(X2)),
             aes(X1, X2, label = cluster))
p_paper

# 
# save(list = c("cl_louvain", "tmp_clusters", "nn_cells", "nn_inothercluster"),
#      file = file.path("~", "asd_analysis", "savepoint", "clusters.RData"))






# Doublets and ambiguous cells ----------------------------------

# load (or re-execute everything in this section):
load(file.path("~", "asd_analysis", "savepoint", "doublets.RData"))







# number of NN from different cluster:
nn_inothercluster <- colSums(
  matrix(tmp_clusters[ t(nn_cells) ],
         ncol = nrow(nn_cells))   != 
  matrix(rep(tmp_clusters, each = ncol(nn_cells)),
         ncol = nrow(nn_cells)) )
 

# in silico doublets: randomly draw cells from different clusters and pool their UMIs to form a "synthetic" doublet:

cellsA <- sample(1:ncol(counts), 50000)
cellsB <- rep(NA, 50000)
smpA <- cellinfo$sample[cellsA]
clA  <- tmp_clusters[cellsA]
tmp <- data.frame(smpA, clA) %>% group_by(smpA, clA) %>% tally
for(i in 1:nrow(tmp)) {
  is_smp <- cellinfo$sample[cellsA] == tmp$smpA[i]  
  is_cl  <- tmp_clusters[cellsA] == tmp$clA[i]
  # sample amongst cells from same sample and different cluster:
  cellsB[ is_smp & is_cl ] <- base::sample(
    x = which(cellinfo$sample == tmp$smpA[i]   & !tmp_clusters == tmp$clA[i]),
    size = tmp$n[i],
    replace = T) # in case one cluster is larger than all others combined
}

doublet_raw <- Ccounts[, cellsA] + Ccounts[, cellsB]
doublet_pcs <- predict(pca,
                       newdata = sqrt( (t(doublet_raw) / colSums(doublet_raw))[, is_informative] ))


# Alternative 1 (clearer):
a <- FNN::get.knn(rbind(pca$x, doublet_pcs), k = 50)
nn_doublets <- a$nn.index
nndists_doublets <- a$nn.dist

# Alternative 2 (faster):
library( RcppAnnoy )
featureMatrix <- rbind(pca$x, doublet_pcs); k_nn <- 50
annoy <- new( AnnoyEuclidean, ncol(featureMatrix) )
for( i in 1:nrow(featureMatrix) ) 
  annoy$addItem( i-1, featureMatrix[i,] )
annoy$build( 50 ) # builds a forest  of n_trees trees. More trees gives higher precision when querying.
nn_doublets <- t( sapply( 1:annoy$getNItems(), function(i) annoy$getNNsByItem( i-1, k_nn) + 1 ) )
nndists_doublets <- sapply( 1:ncol(nn_doublets), function(j) sqrt( rowSums( ( featureMatrix - featureMatrix[ nn_doublets[,j], ] )^2 ) ) )
rm(featureMatrix, annoy)


# percentage of synthetic doublets in neighborhood for each cell:
dblts_perc <- rowMeans( nn_doublets > ncol(counts) )[ 1:ncol(counts) ]



# Run UMAP with Annoy's output
ump2 <- uwot::umap( NULL, nn_method = list( idx=nn_doublets, dist=nndists_doublets), 
                    n_threads=40, spread = 15, verbose=TRUE )

is_synth <- 1:nrow(ump2) > nrow(pca$x)






# save(list = c("nn_doublets", "nndists_doublets", "cellsA", "cellsB",
#                 "dblts_perc", "is_synth", "ump2"),
#      file = file.path("~", "asd_analysis", "savepoint", "doublets.RData"))








  









# DESeq -------------------------------------------------------------------
library(DESeq2)
library(BiocParallel)

# visualize dirty cells we clean away:
tmp <- data.frame(umap_euc,
                  diagnosis = cellinfo$diagnosis,
                  clean = dblts_perc < 3/50  & nn_inothercluster < 1,
                  Gene = Tcounts[, "TTF2"] / sfs/mean(1/sfs),
                  cl = factor(tmp_clusters))
ggplot() + coord_fixed()+
  geom_point(data=filter(tmp,  clean), aes(X1, X2, col = cl), size=.1) +
  geom_point(data=filter(tmp, !clean), aes(X1, X2), col = "black", size=.1) +
  geom_label(data=group_by(tmp, cl) %>% summarise(X1=mean(X1), X2=mean(X2)), aes(X1, X2, label=cl))

tmp <- as.matrix(table(sample=cellinfo$sample, clean = dblts_perc < 3/50  & nn_inothercluster < 1))
data.frame(sample = rownames(tmp), dirtyProportion = tmp[,1] / (tmp[,1] + tmp[,2])) %>% left_join(sampleTable, by="sample") %>% ggplot(aes(sample, dirtyProportion, col = diagnosis))+geom_point()



# compute for a single cluster
sel <-Tcounts[, "SYT1"] > 1 & Tcounts[, "CUX2"] > 0  & dblts_perc < 3/50  & nn_inothercluster < 1
pseudobulks <- as.matrix(t( fac2sparse(cellinfo$sample[sel]) %*% t(Ccounts[, sel]) ))
coldat <- filter(sampleTable, sample %in% colnames(pseudobulks)) %>% 
  mutate(individual = factor(individual),
         diagnosis = factor(diagnosis, levels = c("Control", "ASD")),
         region    = factor(region))
rownames(coldat) <- coldat$sample

dds <- DESeqDataSetFromMatrix( pseudobulks,
                               coldat[colnames(pseudobulks), ],
                               design = ~ sex + region + age + diagnosis )
# For cluster 5, I tested that we do not need interactions between sex, region and diagnosis. I used
# DESeq's LTR for this (see mail to Simon at mid-September 2019).
dds <- DESeq(dds, 
             parallel=TRUE, BPPARAM=MulticoreParam(20))
res_df <- results(dds, name = "diagnosis_ASD_vs_Control") %>% as.data.frame() %>% rownames_to_column("Gene")
table(res_df$padj < .1)  
res_df %>% arrange(padj) %>% head(n=20)





