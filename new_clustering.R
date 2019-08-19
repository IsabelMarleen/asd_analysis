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

raw <- readTENxH5File( "~/sds/sd17l002/u/anders/tmp/ASD.h5" )
rawt <- t(raw)

cs <- colSums(raw)

meta <- read_tsv("~/sds/sd17l002/u/anders/tmp/ASD2/meta.txt" )

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

nn <- ump$nn$cosine$idx[,-1]
ump <- ump$embedding

plot( ump, pch=".", col=adjustcolor( as.integer(as.factor(cellTable$cluster)), .1 ), asp=1 )
text( tapply( ump[,1], cellTable$cluster, median ), tapply( ump[,2], cellTable$cluster, median ), 
      levels(as.factor(cellTable$cluster)) )


# Get Louvain clusters

library( igraph )
g <- graph_from_edgelist( as.matrix( map_dfr( 1:10, function(i) data.frame( from=1:nrow(nn), to=nn[,i] ) ) ) )
louv <- cluster_louvain( as.undirected( g, "collapse" ) )

ncl <- max(membership(louv))

plot( ump, pch=".", col = rainbow( ncl, alpha=.1, v=.7 )[ membership(louv) ], asp=1 )
text( tapply( ump[,1], membership(louv), median ), tapply( ump[,2], membership(louv), median ), 1:ncl )


# Border cells: How many of the k nearest neighbors of a cell are not in the same cluster as the cell?
foreignNN <- rowSums( matrix( membership(louv)[ nn ], ncol=ncol(nn) ) != membership(louv) )
table( foreignNN )

# These cells have no foreign neighbors
interior <- foreignNN == 0

plot( ump, pch=".", col = adjustcolor( ifelse( interior, "black", "red" ), .1 ), asp=1 )

# For which cells is the first nearest neighbor another cluster and which is it?
tibble( cell_cluster = membership(louv), nn1_cluster = membership(louv)[nn[,1]] ) %>%
  group_by_all %>% tally %>% spread( nn1_cluster, n, fill=0 ) %>% 
  column_to_rownames("cell_cluster") -> cluster_nbgh

# Show graph of clusters with at least x% of cells having a neighbor in the other cluster
plot( graph_from_adjacency_matrix( cluster_nbgh / sizes(louv) > .01 ) )

comp <- components( graph_from_adjacency_matrix( cluster_nbgh / sizes(louv) > .02 ) )
gcms <- comp$membership[ membership(louv)  ] # grand cluster membership
foreignNN <- rowSums( matrix( gcms[ nn ], ncol=ncol(nn) ) !=gcms )

plot( ump, pch=".", col = ifelse( foreignNN>0, "#00000018", rainbow( ncl, alpha=.1, v=.7 )[ gcms ] ), asp=1 ) 
text( tapply( ump[,1], gcms, median ), tapply( ump[,2], gcms, median ), 1:20 )

#Calculation of pseudobulks for one cluster at a time
#Add new clustering to cellTable as new column

cellTable <- cellTable %>%
  mutate(newcluster = sprintf( "MC%02d", gcms) )

counts <- TENxMatrix( "/home/anders/pub/ASD.h5", "matrix" )
do_DESeq_on_newcluster <- function( cluster ){
  pseudobulk <- sapply( sampleTable$sample, function(s)
    rowSums( counts[ , meta$sample == s & cellTable$newcluster==cluster, drop=FALSE ] ) )
  dds <- DESeqDataSetFromMatrix( pseudobulk, sampleTable, ~ sex+region+age+diagnosis)
  dds <- dds[ rowSums(counts(dds)) >= 10, colSums(counts(dds)) > 0 ]
  dds <- DESeq( dds )
  dds
}

plan( multiprocess, workers=20 )
dds_nc <- future_map( sort( setdiff( unique(cellTable$newcluster), "MC18") ), do_DESeq_on_newcluster )
names(dds_nc) <- sort( setdiff( unique(cellTable$newcluster), "MC18") )
      
genenames <- tibble( ensg = h5read("ASD.h5", "matrix/genes"), 
        name = h5read("ASD.h5", "matrix/gene_names") )
res.paper <- read.delim("~/sds/sd17l002/u/isabel/S4.csv", sep=";", dec = ",") %>%
  select(cluster = Cell.type, gene = gene.ID, name=Gene.name, q.value)

res.tibble <- map_dfr( results, as_tibble, rownames="gene", .id="cluster" ) %>%
  left_join( genenames, by=c( "gene" = "ensg" ) ) #%>%
  #filter( .$padj < .1 )

res.nc <- map_dfr( dds_nc, function(x) x %>% results %>% as_tibble( rownames="gene" ), .id="newcluster" ) %>%
  left_join( genenames, by=c( "gene" = "ensg" ) )
 

#Plot that relates significance from DESeq output to MAST output
res.tibble %>%
  left_join(res.paper)%>%
  mutate(signif.DESeq = !is.na(padj) & padj < .05, 
         signif.MAST = !is.na(q.value) & q.value <.05) %>% 
  ggplot+
      geom_point(aes(signif.MAST, signif.DESeq, col=signif.DESeq), size=.1, position= "jitter")+
      facet_wrap(~cluster)

#Plot that compares DESeq on manual clusters and paper clusters
#Currently this plots all genes, which takes long and is not very smart
map_dfr( results_nc, as_tibble, rownames="gene", .id="newcluster" ) %>%
  mutate(padjnc = .$padj,
            padj=NULL) %>%
  left_join(res.tibble, by="gene")%>%
  mutate(sign.nc = padjnc < .1, sign.pc = padj <.1) %>% 
ggplot+
  geom_point(aes(sign.nc, sign.pc, col=sign.nc), size=.1, position= "jitter")+
  facet_wrap(~newcluster)

#look at diagnostic marker genes for different EN layers
  ggplot()+
    geom_point(aes( ump[, 1], ump[,2] ,col= raw["TLE4", ]/cs), size=.01 )+
    scale_color_gradientn(colours = rev(rje::cubeHelix(100)), trans="sqrt")+
    geom_text( aes( tapply( ump[,1], cellTable$cluster, median ), 
                    tapply( ump[,2], cellTable$cluster, median ), 
        label=levels(factor(cellTable$cluster)) ), data=NULL )
  
  
  
#Plot for Bag3 ASD vs Normal
k <-   tibble(NONO = raw["NONO", ], cell=meta$cell, fracMT=raw["MTND2P28", ]/cs , 
              BAG3 = raw["BAG3", ], fracNONO= raw["NONO", ]/cs, newcluster=cellTable$newcluster)%>%
         left_join(meta, by="cell")

  
ggplot(k, aes(k$fracMT, k$diagnosis, col=k$diagnosis))+
  geom_point( position = "jitter", size=.1)+
  facet_wrap(~k$sample)
  
  ggplot(k, aes(sqrt(k$fracNONO), k$sample, col=k$diagnosis))+
    geom_point( position = "jitter", size=.1)+
    facet_wrap(~k$newcluster)
  
  k %>% filter(newcluster=="MC08") %>% group_by(sample, diagnosis) %>% 
    summarise(y=mean(sqrt(fracNONO)>0.01) ) %>% ggplot() + geom_point(aes(x=diagnosis, y=y ))
