
sampleTable <-
  cellinfo %>% select( sample : RNA.Integrity.Number ) %>% unique
sampleTable

pseudobulk_L23 <-
  sapply( sampleTable$sample, function(s)
    rowSums( counts[ , cellinfo$sample == s & cellinfo$cluster ==  drop=FALSE ] ) )
pseudobulk_L23[ 1:5, 1:5 ]


formula1 <- ~   Seqbatch + 
  region + age + sex +  RNA.Integrity.Number + post.mortem.interval..hours. + diagnosis
formula2 <- ~  Seqbatch + 
  region + age + sex +  RNA.Integrity.Number + diagnosis
formula3 <- ~ Seqbatch + 
  region + age + sex + diagnosis
  
formula4 <- ~diagnosis

do_DESeq <- function( formula ){
  
  dds2 <- DESeqDataSetFromMatrix( pseudobulk_L23, sampleTable, formula )
  dds2 <- estimateSizeFactors(dds2)
  nc <- counts(dds2, normalized=TRUE)
  filter <- rowSums(nc >= 10) >= 2
  dds2 <- dds2[filter,]
  #dds2 <- DESeq( dds2 )
  dds2 <- estimateDispersions(dds2)
  dds2 <- nbinomWaldTest(dds2, maxit=500)
  return(results( dds2 ))
}

do_DESeq2 <- function(formula) {
  dds <- DESeqDataSetFromMatrix( pseudobulk_L23, sampleTable, formula )
  dds <- DESeq( dds )
  return(results( dds ))
}




min <- do_DESeq(formula4)
min2 <- do_DESeq2(formula4)
full <-  do_DESeq(formula1)
oneoff_PMI <- do_DESeq(formula2)
twooff_PMI_RIN <- do_DESeq(formula3)
