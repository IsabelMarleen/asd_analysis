# for pre-processing scRNAseq data:
library(Matrix)
library(irlba)
library(uwot)
library(FNN)
library(RcppAnnoy)
library(igraph)
# 
library(tidyverse)
library(cowplot)
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

path <- "~/sds/sd17l002/p/MS/"


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
            select(cell, tSNE1=tsne1, tSNE2=tsne2, everything())
# info per patient:
sampleTable <-
  cellinfo %>% select( sample : RIN, lesion_stage ) %>% unique
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
set.seed(100) # seed ensures that UMAP gives reproducible result
pca <- irlba::prcomp_irlba( x = sqrt(t(norm_counts[is_informative,])),
                            n = 40,
                            scale. = TRUE)

# Find nearest neighbors for UMAP and Louvain clustering:
set.seed(100) # seed ensures that UMAP gives reproducible result
featureMatrix <- pca$x; k_nn <- 50
annoy <- new( AnnoyEuclidean, ncol(featureMatrix) )
for( i in 1:nrow(featureMatrix) ) 
  annoy$addItem( i-1, featureMatrix[i,] )
annoy$build( 50 ) # builds a forest  of n_trees trees. More trees gives higher precision when querying.
nn_cells <- t( sapply( 1:annoy$getNItems(), function(i) annoy$getNNsByItem( i-1, k_nn) + 1 ) )
nndists_cells <- sapply( 1:ncol(nn_cells), function(j) sqrt( rowSums( ( featureMatrix - featureMatrix[ nn_cells[,j], ] )^2 ) ) )
rm(featureMatrix, annoy)


set.seed(100) # for reproducibility, we also have to set n_sgd_threads to 1.
u <- uwot::umap( pca$x, spread = 10, n_threads = 1, n_sgd_threads=1,
                 nn_method = list(idx=nn_cells, dist=nndists_cells)) 
u <- as_tibble(u, .name_repair = ~ c("u1", "u2"))


# simplified paper_clusters in UMAP
paper_clusters <- case_when(
  grepl("^EN-", cellinfo$cell_type)  ~ "excit_layers", 
  grepl("^IN-", cellinfo$cell_type)  ~ "inhibitory", 
  grepl("^OL-", cellinfo$cell_type)  ~ "OL", 
  grepl("^AST", cellinfo$cell_type)  ~ "Astro", 
  TRUE ~ cellinfo$cell_type) %>% factor()
ggplot() +
  geom_point(data = u %>% bind_cols(paper_clusters = paper_clusters),
             aes(u1, u2, col = paper_clusters), size=.5) + coord_fixed() +
  geom_label(data = u %>% bind_cols(paper_clusters = paper_clusters) %>%
               group_by(paper_clusters) %>% summarise(u1=mean(u1), u2=mean(u2)),
             aes(u1, u2, label = paper_clusters)) + theme(legend.position = "none") 



# Louvain Clustering ------------------------------------------------------

# We use the approximate nearest neighbors computed above.

# matrix has to be sparse, otherwise takes 80 GB of RAM:
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
  # Returns normalized values of the form k/s/mean(1/s), where
  # k are UMIs and s are a cell's total UMIs across all genes (aka library size).
  # Side effect: saves smoothing to knn_list so we don't have to recompute.
  if(is.null(knn_list[[g]])){ 
    norm_umis <- matrix(Ccounts[g, c(nn_cells)] / sfs[c(nn_cells)],
                        ncol = ncol(nn_cells))
    knn_smoothed <- rowSums(norm_umis * t(kernel_weights)) / colSums(kernel_weights)
    knn_list[[g]] <<- knn_smoothed / mean(1/sfs)
  }
  
  knn_list[[g]]
}





# Marker examples ---------------------------------------------------------



major_celltypes <- matrix(0, ncol=9, nrow=8); diag(major_celltypes) <- 1
dimnames(major_celltypes) <- list(Class = c("Oligod", "OPC", "Astrocyte",
                                              "Microglia", "endo",
                                              "nNeuron","iNeuron", "eNeuron" ),
                                    Gene  = c("PLP1", "TNR", "AQP4", 
                                              "PTPRC_ENSG00000081237",
                                              "VWF", "SYT1", "GAD1", "SATB2", "NRGN"))
major_celltypes[c("eNeuron","iNeuron"), "SYT1"] <- 1
major_celltypes["nNeuron", "NRGN"] <- 1





# Classification ----------------------------------------------------------

# plotting function:
em_result_plot <- function(probs=p, ump =u[sel,], p_thresh = .5) {
  p_umap <- ump %>%
    bind_cols(class=apply(probs, 1, function(x)
      ifelse(max(x) > p_thresh, colnames(probs)[which.max(x)], NA) )) %>%
    mutate(class = factor(class, levels = colnames(probs))) %>%
    ggplot(aes(u1, u2, col = class)) + geom_point(size=.5) + coord_fixed() +
    scale_color_manual(values = scales::hue_pal()(ncol(probs)),
                       drop=F, na.value = "grey") 
  return(p_umap)
}

# smooth and visualize (takes a while on first computation):
data.frame(sapply(colnames(major_celltypes), knn_smooth)) %>%
  bind_cols(u) %>%
  gather(Gene, knn, -u1, -u2) %>%
  ggplot(aes(u1, u2, col = knn)) + geom_point(size=.5)+coord_fixed()+
  facet_wrap(~Gene) + scale_color_sqrt()


# set-up testbed:
ms <- major_celltypes
expr <- data.frame(sapply(colnames(ms),
                          knn_smooth)  + .1/50)  # +.1/50 avoids zeros
sel <- TRUE # take ALL cells, i.e. 104k
# sel <- male_controls_pfc
expr <- expr[sel,]


p <-learnClasses(ms, expr)
p_bak <- p
em_result_plot(p, p_thresh = .5)







# manual EM ---------------------------------------------------------------

marker_table <- major_celltypes
expr_table <- expr

# number of iterations:
iter = 0
maxiter = 1000
# other housekeeping:
loglik = rep(-Inf, nrow(marker_table)) 
delta = +Inf
tolerance = 0.001
miniter = 10

# initialize p
p <- scpr:::priors_geometric(marker_table, expr_table)
logp <- log(p)

while ((delta > tolerance) && (iter <= maxiter) || (iter < miniter)) {
  
  p <- pmax(exp(logp), 1e-05)
  p <- p/(rowSums(p)+1)
  logp <- sapply(1:nrow(marker_table), function(class) {
    loglik_mat <- sapply(1:ncol(marker_table), function(gene) {
      probs <- p[, class]
      return(scpr:::Lgamma(expr_table[, gene], probs, 
                           log = T))
    })
    log(mean(p[, class])) + rowSums(loglik_mat)
  })
  logp <- logp - log(rowSums(exp(logp))+1)
  loglikOld = loglik
  loglik <- colSums(logp)
  delta <- max(abs(loglikOld - loglik))
  iter = iter + 1
}
p <- exp(logp)
colnames(p) <- rownames(marker_table)


# seems like we are excluding doublets. Still, 90 % of
# excluded cells come from MS, which is clearly artificial and
# needs more work:
classes <- apply(p1, 1, function(x) 
  ifelse(max(x) > .5, colnames(p1)[which.max(x)], NA) )
table(na=is.na(classes), treat = cellinfo$diagnosis)



# -ncol(p) excludes "other" class
scale_hist <- function(x = expr[, "SYT1"], probs=p[, -ncol(p)], highlight = "eNeuron"){
ggplot() +
  geom_histogram(data = data.frame(Gene = x),
                 aes(Gene, stat(density)), bins=100)+
  scale_x_log10(limits = c(min(x), max(x))) +
  # add all classes (dashed):
  lapply(as.data.frame(probs), function(ps){
    param <- scpr::mom_gamma(x, ps)
    geom_density(data=data.frame(theo = rgamma(10000, param[1], param[2])),
                 aes(theo, mean(ps)*stat(density)), linetype = "dashed")})  +
  # emphasize class in red:
  geom_density(
    data=data.frame(theo = rgamma(10000,
                                  scpr::mom_gamma(x, probs[, highlight])[1],
                                  scpr::mom_gamma(x, probs[, highlight])[2])),
    mapping = aes(theo, mean(probs[,highlight])*stat(density)), col = "red") 
}



scale_hist(expr[,"SYT1"], highlight = "eNeuron")













