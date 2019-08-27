#Setup
library( ggplot2 )
library( Matrix )
library( tidyverse )
source( "spmvar.R" )
library( HDF5Array )
library(DESeq2)
library( furrr )

#Load data
path <- "~/sds/sd17l002/u/isabel/rawMatrix/"
cellinfo <- read.delim( paste0( path, "meta.txt" ), stringsAsFactors=FALSE )

# sample table
sampleTable <-
  cellinfo %>% select( sample : RNA.Integrity.Number ) %>% unique

sampleTable$diagnosis <- relevel( factor( sampleTable$diagnosis ), "Control" )

counts <- TENxMatrix( "/home/anders/pub/ASD.h5", "matrix" )
do_DESeq_on_cluster <- function( cluster ){
  pseudobulk <- sapply( sampleTable$sample, function(s)
    rowSums( counts[ , meta$sample == s & cellTable$cluster==cluster, drop=FALSE ] ) )
  dds <- DESeqDataSetFromMatrix( pseudobulk, sampleTable, ~ age+region+sex+diagnosis )
  dds <- dds[ rowSums(counts(dds)) >= 10, colSums( counts( dds ) ) > 0 ]
  dds <- DESeq( dds )
  results( dds )
}

plan( multiprocess, workers=20 )
dds <- future_map( unique( cellTable$cluster ), do_DESeq_on_cluster )
names( dds ) <- unique( cellTable$cluster )

map_dbl( results, ~ sum( .$padj<.1, na.rm=TRUE ) )

res.tibble %>%
  group_by( gene ) %>%
  summarise( nclsig = sum( padj < .1, na.rm=TRUE ) ) %>%
  arrange( -nclsig )

#How many genes below FDR of 10% in each cluster?
#Paper clusters
res.tibble %>%
  group_by( cluster ) %>%
  summarise( ngincl = sum( padj < .1, na.rm=TRUE ) ) %>% #no. of signif genes in cluster
  arrange( -ngincl )

#Manual Clusters
res.nc %>%
  left_join(clusterkey) %>%
  group_by( newcluster, clusterkey ) %>%
  summarise( nginncl = sum( padj < .1, na.rm=TRUE ) ) %>% #no. of signif. genes in new cluster
  arrange( -nginncl )

#Table of padj of genes in DESeq and MAST
tibble( ensg = h5read("ASD.h5", "matrix/genes"), 
        name = h5read("ASD.h5", "matrix/gene_names") )

