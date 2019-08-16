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
  dds <- DESeqDataSetFromMatrix( pseudobulk, sampleTable, ~ diagnosis)
  dds <- dds[ rowSums(counts(dds)) >= 10, colSums(counts(dds)) > 0 ]
  dds <- DESeq( dds )
  results(dds)
}

plan( multiprocess, workers=20 )
results <- future_map( unique(cellTable$cluster), do_DESeq_on_cluster )
names(results) <- unique(cellTable$cluster)

map_dbl( results, ~ sum( .$padj<.1, na.rm=TRUE ) )

map_dfr( results, as_tibble, rownames="gene", .id="cluster" ) %>%
  group_by( gene ) %>%
  summarise( nclsig = sum( padj < .1, na.rm=TRUE ) ) %>%
  arrange( -nclsig )

#How many genes below FDR of 10% in each cluster?
#Paper clusters
map_dfr( results, as_tibble, rownames="gene", .id="cluster" ) %>%
  group_by( cluster ) %>%
  summarise( ngincl = sum( padj < .1, na.rm=TRUE ) ) %>% #no. of genes in cluster
  arrange( -ngincl )

#Manual Clusters
map_dfr( results_nc, as_tibble, rownames="gene", .id="newcluster" ) %>%
  group_by( newcluster ) %>%
  summarise( nginncl = sum( padj < .1, na.rm=TRUE ) ) %>% #no. of genes in new cluster
  arrange( -nginncl )

#Table of padj of genes in DESeq and MAST
tibble( ensg = h5read("ASD.h5", "matrix/genes"), 
        name = h5read("ASD.h5", "matrix/gene_names") )

