# for pre-processing scRNAseq data:
library(Matrix)
library(irlba)
library(uwot)
library(FNN)
library(RcppAnnoy)
library(igraph)
# 
library(tidyverse)
library(scpr)
library(ggpower) # for power-transformations in ggplot
library(rje)



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



# read in data ------------------------------------------------------------
# The matrix with raw UMI counts can be downloaded from [here](http://cells.ucsc.edu/?ds=autism#).
# This section reads them into R.
                                                             
path <- "~/sds/sd17l002/p/ASD/"


counts <- readMM( file.path( path, "rawMatrix", "matrix.mtx" ) )
# make gene symbols unique (by concatenating ensembleID where necessary):
gene_info <- read.delim( file.path( path, "rawMatrix", "genes.tsv" ), header=FALSE, as.is=TRUE ) %>%
  mutate(unique = case_when(
    duplicated(V2) | duplicated(V2, fromLast=T) ~ paste(V2, V1, sep="_"),
    TRUE ~ V2))
rownames(counts) <- gene_info$unique
colnames(counts) <- readLines( file.path( path, "rawMatrix", "barcodes.tsv" ) )

# info per cell (append tSNE we downloaded from cells.ucsc.edu as well):
cellinfo <- read.delim( file.path( path, "rawMatrix", "meta.txt" ), stringsAsFactors=FALSE ) %>%
  left_join(read_tsv(file.path(path, "rawMatrix", "tSNE.coords.tsv"),
                     col_names = c("cell", "tSNE1", "tSNE2")),
            by = "cell") %>% select(cell, tSNE1, tSNE2, everything())
# info per patient:
sampleTable <-
  cellinfo %>% select( sample : RNA.Integrity.Number ) %>% unique
sampleTable



# PCA and UMAP ------------------------------------------------------------

# if RAM/memory is not an issue, this speeds up downstream analysis:
Tcounts <- as(t(counts), "dgCMatrix") #  fast: Tcounts[, "SYN1"]
Ccounts <- as(counts, "dgCMatrix")    #  fast: Ccounts[, 1337] & colSums(Ccounts)


sfs <- colSums(Ccounts)
norm_counts <- t(t(Ccounts) / sfs)
rownames(norm_counts) <- rownames(Ccounts)
poisson_vmr <- mean(1/sfs)
gene_means <- rowMeans( norm_counts )
gene_vars <- rowVars_spm( norm_counts )
is_expressed <- colSums( Tcounts != 0 ) > 100
is_informative <- gene_vars/gene_means > 1.5 * poisson_vmr  &  is_expressed
plot(gene_means, gene_vars/gene_means, pch=".", log = "xy")
points(gene_means[is_informative], (gene_vars/gene_means)[is_informative], pch=".", col = "red" )

pca <- irlba::prcomp_irlba( x = sqrt(t(norm_counts[is_informative,])),
                            n = 40,
                            scale. = TRUE)
seed(100)
u <- uwot::umap( pca$x, spread = 10, n_threads = 40) # euc: euclidean distance
u <- as_tibble(u, .name_repair = ~ c("u1", "u2"))



# Louvain Clustering ------------------------------------------------------
set.seed(100)
featureMatrix <- pca$x; k_nn <- 50
annoy <- new( AnnoyEuclidean, ncol(featureMatrix) )
for( i in 1:nrow(featureMatrix) ) 
  annoy$addItem( i-1, featureMatrix[i,] )
annoy$build( 50 ) # builds a forest  of n_trees trees. More trees gives higher precision when querying.
nn_cells <- t( sapply( 1:annoy$getNItems(), function(i) annoy$getNNsByItem( i-1, k_nn) + 1 ) )
nndists_cells <- sapply( 1:ncol(nn_cells), function(j) sqrt( rowSums( ( featureMatrix - featureMatrix[ nn_cells[,j], ] )^2 ) ) )
rm(featureMatrix, annoy)


# has to be sparse, otherwise takes 80 GB of RAM:
adj <- Matrix(0, nrow = nrow(pca$x), ncol = nrow(pca$x)) 
for(i in 1:ncol(nn_cells))
  adj[ cbind(1:nrow(pca$x), nn_cells[, i]) ] <- 1
for(i in 1:ncol(nn_cells))
  adj[ cbind(nn_cells[, i], 1:nrow(pca$x)) ] <- 1
set.seed(100)
cl_louvain <- cluster_louvain(  graph_from_adjacency_matrix(adj, mode = "undirected") )
cl_louvain <- factor(cl_louvain$membership)



# set up smoothing --------------------------------------------


scpr_list <- list()
scpr_smooth <- function(g= "SYN1"){
  # Smooth gene expression with scpr.
  # Side effect: saves smoothing to scpr_list so we don't have to recompute.
  if(is.null(scpr_list[[g]])){
    scpr_list[[g]] <<- scpr(Tcounts[, g], sfs, pca$x, lambda=1,
                            nn_list = list(nn_cells, nndists_cells)) / sfs/mean(1/sfs)
  }
  scpr_list[[g]]
}





kernel_weights <- apply(nndists_cells,
                        1,
                        function(d) scpr::tricube_fromNN1(d, max(d)))

knn_list <- list()
knn_smooth <- function(g = "SYT1"){
  # knn-smooth a gene with tricube weights.
  # Side effect: saves smoothing to knn_list so we don't have to recompute.
  if(is.null(knn_list[[g]])){ 
    norm_umis <- matrix(Ccounts[g, c(nn_cells)] / sfs[c(nn_cells)],
                        ncol = ncol(nn_cells))
    knn_smoothed <- rowSums(norm_umis * t(kernel_weights)) / colSums(kernel_weights)
    knn_list[[g]] <<- knn_smoothed / mean(1/sfs)
  }
  
  knn_list[[g]]
}





# Classification ----------------------------------------------------------



sapply(c(
"SLC17A7", # excitatory neurons
"PLP1",    # oligodendrocytes
"SLC1A2"  # glial / astro
), knn_smooth)


data.frame(sapply(c(
  "SLC17A7", # excitatory neurons
  "PLP1",    # oligodendrocytes
  "SLC1A2"  # glial / astro
), knn_smooth)) %>% bind_cols(u) %>%
  gather(Gene, knn, -u1, -u2) %>%
  ggplot(aes(u1, u2, col = knn)) + geom_point(size=.5)+coord_fixed()+
  facet_wrap(~Gene) + scale_color_sqrt()


