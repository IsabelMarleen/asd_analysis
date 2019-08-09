#Setup
library( ggplot2 )
library( Matrix )
library( tidyverse )
source( "spmvar.R" )
library( HDF5Array )

library( furrr )

#Load data
path <- "~/sds/sd17l002/u/isabel/rawMatrix/"
cellinfo <- read.delim( paste0( path, "meta.txt" ), stringsAsFactors=FALSE )

# sample table
sampleTable <-
  cellinfo %>% select( sample : RNA.Integrity.Number ) %>% unique

sampleTable$diagnosis <- relevel( factor( sampleTable$diagnosis ), "Control" )

cluster <- "Microglia"
do_DESeq_on_cluster <- function( cluster ){
  counts <- TENxMatrix( "~/sds/sd17l002/u/anders/tmp/ASD.h5", "matrix" )
  pseudobulk <- sapply( sampleTable$sample, function(s)
    rowSums( counts[ , cellinfo$sample == s & cellinfo$cluster==cluster, drop=FALSE ] ) )
  dds <- DESeqDataSetFromMatrix( pseudobulk, sampleTable, ~ diagnosis)
  dds <- dds[ rowSums(counts(dds)) >= 10, colSums(counts(dds)) > 0 ]
  dds <- DESeq( dds )
  results(dds)
}

plan( multiprocess, workers=20 )
results <- future_map( unique(cellinfo$cluster), do_DESeq_on_cluster )
names(results) <- unique(cellinfo$cluster)

map_dbl( results, ~ sum( .$padj<.1, na.rm=TRUE ) )

map_dfr( results, as_tibble, rownames="gene", .id="cluster" ) %>%
  group_by( gene ) %>%
  summarise( nclsig = sum( padj < .1, na.rm=TRUE ) ) %>%
  arrange( -nclsig )

?as_tibble
