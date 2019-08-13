library( Matrix )
library( tidyverse )
library( rhdf5 )
library( HDF5Array )

readTENxH5File <- function( filename )
  sparseMatrix(
    p = h5read( filename, "matrix/indptr" ),
    i = h5read( filename, "matrix/indices" ) + 1,
    x = as.numeric( h5read( filename, "matrix/data" ) ),
    dims = h5read( filename, "matrix/shape" ),
    dimnames = list(
      as.vector( h5read( filename, "matrix/gene_names" ) ),
      as.vector( h5read( filename, "matrix/barcodes" ) )
    )
  )

raw <- readTENxH5File( "~/tmp/ASD/ASD.h5" )
rawt <- t(raw)

cs <- colSums(raw)

meta <- read_tsv("~/tmp/ASD/meta.txt" )

meta %>% 
select( sample : `RNA Integrity Number` ) %>%
unique ->
  sampleTable

meta %>% 
select( -( individual : `RNA Integrity Number` ) )  ->
  cellTable

stopifnot( all( str_sub( cellTable$cell, 1, 14 ) == colnames(raw) ) )
stopifnot( all( str_sub( cellTable$cell, 1, 14 ) == rownames(rawt) ) )

colVars_spm <- function( spm ) {
  stopifnot( is( spm, "dgCMatrix" ) )
  ans <- sapply( seq.int(spm@Dim[2]), function(j) {
    mean <- sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
    sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) +
      mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) } ) / ( spm@Dim[1] - 1 )
  names(ans) <- spm@Dimnames[[2]]
  ans
}

gene_means <- colMeans( rawt/cs )
gene_vars <- colVars_spm( rawt/cs )
xi <- mean( 1/cs )

plot( gene_means, gene_vars / gene_means, pch=".", log="xy" )
abline( h = 1.5*xi )

infrmgenes <- gene_vars / gene_means > 1.5*xi & gene_means > 3e-8
rawi <- raw[ infrmgenes, ]

pca <- irlba::prcomp_irlba( sqrt(t(rawi)/cs), 30 )
ump <- uwot::umap( pca$x, metric="cosine", ret_nn=TRUE, n_threads=40, min_dist=0.03, verbose=TRUE )

nn <- ump$nn$cosine[,-1]
ump <- ump$embedding

plot( ump, pch=".", col=adjustcolor( as.integer(as.factor(cellTable$cluster)), .1 ), asp=1 )
text( tapply( ump[,1], cellTable$cluster, median ), tapply( ump[,2], cellTable$cluster, median ), levels(as.factor(cellTable$cluster)) )


# Get Louvain clustes

library( igraph )
g <- graph_from_edgelist( as.matrix( map_dfr( 1:10, function(i) data.frame( from=1:nrow(nn$idx), to=nn$idx[,i] ) ) ) )
louv <- cluster_louvain( as.undirected( g, "collapse" ) )

ncl <- max(membership(louv))

plot( ump, pch=".", col = rainbow( ncl, alpha=.1, v=.7 )[ membership(louv) ], asp=1 )
text( tapply( ump[,1], membership(louv), median ), tapply( ump[,2], membership(louv), median ), 1:ncl )


# Border cells: How many of the k nearest neighbors of a cell are not in the same cluster as the cell?
foreignNN <- rowSums( matrix( membership(louv)[ nn$idx ], ncol=ncol(nn$idx) ) != membership(louv) )
table( foreignNN )

# These cells have no foreign neighbors
interior <- foreignNN == 0

plot( ump, pch=".", col = adjustcolor( ifelse( interior, "black", "red" ), .1 ), asp=1 )

# For which cells is the first nearest neighbor another cluster and which is it?
tibble( cell_cluster = membership(louv), nn1_cluster = membership(louv)[nn$idx[,1]] ) %>%
  group_by_all %>% tally %>% spread( nn1_cluster, n, fill=0 ) %>% column_to_rownames("cell_cluster") -> cluster_nbgh

# Show graph of clusters with at least x% of cells having a neighbor in the other cluster
plot( graph_from_adjacency_matrix( cluster_nbgh / sizes(louv) > .01 ) )

comp <- components( graph_from_adjacency_matrix( cluster_nbgh / sizes(louv) > .02 ) )
gcms <- comp$membership[ membership(louv)  ] # grand cluster membership
foreignNN <- rowSums( matrix( gcms[ nn$idx ], ncol=ncol(nn$idx) ) !=gcms )

plot( ump, pch=".", col = ifelse( foreignNN>0, "#00000018", rainbow( ncl, alpha=.1, v=.7 )[ gcms ] ), asp=1 ) 

