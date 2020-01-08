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

# Find nearest neighbors for UMAP and Louvain clustering:
set.seed(100)
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
  grepl("^L", cellinfo$cluster)  ~ "excit_layers", 
  grepl("^IN", cellinfo$cluster)  ~ "inhibitory", 
  grepl("^AST", cellinfo$cluster)  ~ "Astro", 
  grepl("^Neu-NRGN", cellinfo$cluster)  ~ "NRGN", 
  TRUE ~ cellinfo$cluster) %>% factor()
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


# threeCelltypes 
three_celltypes_1 <- matrix(0, ncol=3, nrow=4); diag(three_celltypes_1) <- 1
dimnames(three_celltypes_1) <- list(Class = c("Oligod", "OPC", "Astrocyte", "other"),
                                  Gene  = c("CNP", "PDGFRA", "AQP4"))
# threeCelltypes - alternative 2:
three_celltypes_2 <- matrix(0, ncol=3, nrow=4); diag(three_celltypes_2) <- 1
dimnames(three_celltypes_2) <- list(Class = c("Oligod", "OPC", "Astrocyte", "other"),
                                  Gene  = c("PLP1", "TNR", "AQP4"))

# inhibitory_exhibitory Neurons:
ie_neurons <- matrix(c(1,1,0, 1,0,0, 0,1,0), nrow=3, dimnames = list(
  Class = c("inhibitory", "excitatory", "other"),
  Gene  = c("SYT1", "GAD1", "SATB2")
))


# fiveCelltypes:
five_celltypes <- matrix(0, nrow=6, ncol = 5); diag(five_celltypes) <- 1
dimnames(five_celltypes) <- list(
  Class = c("Oligod", "OPC", "Astro", "Microglia", "endo", "other"),
  Gene  = c("PLP1", "TNR", "AQP4", "PTPRC_ENSG00000081237", "VWF")
)

five_plusN <- matrix(0, nrow=7, ncol = 6); diag(five_plusN) <- 1
dimnames(five_plusN) <- list(
  Class = c("Oligod", "OPC", "Astro", "Microglia", "endo", "Neuron", "other"),
  Gene  = c("PLP1", "TNR", "AQP4", "PTPRC_ENSG00000081237", "VWF", "SYT1")
)

six_noOther <- five_plusN; six_noOther <- six_noOther[-7,]


five_plusN_subtypes <- matrix(0, nrow=8, ncol = 8); diag(five_plusN_subtypes) <- 1
dimnames(five_plusN_subtypes) <- list(
  Class = c("Oligod", "OPC", "Astro", "Microglia", "endo", "iNeuron", "eNeuron", "other"),
  Gene  = c("PLP1", "TNR", "AQP4", "PTPRC_ENSG00000081237", "VWF",
            "SYT1", "GAD1", "SATB2")
)
five_plusN_subtypes["iNeuron", "GAD1"] <- 1
five_plusN_subtypes["eNeuron", "GAD1"] <- 0
five_plusN_subtypes["eNeuron", c("SYT1","SATB2")] <- 1
five_plusN_subtypes["other", "SATB2"] <- 0

# leave out other:
five_plusN_subtypes_noOther <- five_plusN_subtypes[-8,]


# specify NRGN Neurons as well:
five_and_threeNeurons <- cbind(five_plusN_subtypes_noOther, NRGN=0)
five_and_threeNeurons <- rbind(five_and_threeNeurons,
                               "nNeurons"=c(0,0,0,0,0,1,0,0,1))


# Literature search (cumbersome) ------------------------------------------
# Markers come from these sources:
#  Velmeshev, Science 2019
#  Amel from Schirmer group ("Marker genes" email from August 2019)
#  Mariani, ..., Vaccarino (2012 in PNAS)
#  Kang, Nature 2011
#    Note that Kang et al. have NOT sorted celltypes, they just curated literature
#    for markers.


# Markers from Kang and Mariani (deleted unexpressed markers):
genes_l5 <- c("BCL11B", "ETV1", "TLE4")])  # note: BCL11B aka CTIP2
# Markers from Kang: (deleted unexpressed markers):
genes_l6 <- c("ZFPM2", "TLE4", "FOXP2", "TBR1", "SOX5")
genes_l234 <- c("CUX1", "UNC5D")
genes_l4   <- c("RORB")

# AQP4 is from Amel, rest from Kang:
genes_astro <- c("GFAP", "S100B", "ALDOC", "AQP4")

#  PTPRC from Amel, rest from Kang:
genes_microglia <- c("CFH", "FCER1G", "TNIP2", "PTPRC_ENSG00000081237")


# Classification ----------------------------------------------------------

# plotting function:
em_result_plot <- function(probs=p, ump =u[sel,], p_thresh = .5) {
 p_umap <- ump %>%
  bind_cols(class=apply(probs, 1, function(x)
    ifelse(max(x) > p_thresh, colnames(probs)[which.max(x)], NA) )) %>%
  mutate(class = factor(class, levels = colnames(probs))) %>%
  ggplot(aes(u1, u2, col = class)) + geom_point(size=.5) +
  scale_color_manual(values = scales::hue_pal()(ncol(probs)),
                     drop=F, na.value = "grey") 
 return(p_umap)
}

# smooth and visualize (takes a while on first computation):
data.frame(sapply(colnames(five_plusN_subtypes), knn_smooth)) %>%
  bind_cols(u) %>%
  gather(Gene, knn, -u1, -u2) %>%
  ggplot(aes(u1, u2, col = knn)) + geom_point(size=.5)+coord_fixed()+
  facet_wrap(~Gene) + scale_color_sqrt()


# some PFC samples with good number of cells:   (PFC was arbitrarily chosen over ACC)
a_male_control <- cellinfo$sample == "5958_BA9"
b_male_asd     <- cellinfo$sample == "5864_BA9"
male_controls_pfc <- cellinfo$region=="PFC" &
  cellinfo$sex == "M" &
  cellinfo$diagnosis == "Control"
male_asd_pfc <- cellinfo$region=="PFC" &
  cellinfo$sex == "M" &
  cellinfo$diagnosis == "ASD"
female_controls_pfc <- cellinfo$region=="PFC" &
  cellinfo$sex == "F" &
  cellinfo$diagnosis == "Control"


# set-up testbed:
# ms <- five_plusN
ms <- five_and_threeNeurons
expr <- data.frame(sapply(colnames(ms),
                          knn_smooth)  + .1/50)
sel <- TRUE # take ALL cells, i.e. 104k
# sel <- male_controls_pfc
expr <- expr[sel,]


p <-learnClasses(ms, expr)
p_bak <- p
em_result_plot(p, p_thresh = .5)














# Insights ----------------------------------------------------------------



# I need coord_trans, not scale_x_* - see ggpower-vignette.

# three_celltypes_2: first EM-iteration fits PLP1 such that Oligods
# become "other" in the next steps.
# Is this a practical problem of shape/rate estimation? Or
# do some "other" cells express PLP1 highly, then this would be
# a theoretical problem of having a marker-less class.

# what's going wrong in five_celltypes? Look again now that you've
# corrected em_histograms function (uses coord_trans now).




