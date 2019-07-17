#Workflow zlm

#Setup
library( ggplot2 )
library( Matrix )
library( locfit )
library( purrr )
library( furrr )
library( tidyverse )
library(MAST)
library(lme4)
source( "spmvar.R" )

#Load data
cellinfo <- read.delim( "~/sds/sd17l002/u/isabel/rawMatrix/meta.txt", stringsAsFactors=FALSE )
counts <- Seurat::Read10X("~/sds/sd17l002/u/isabel/rawMatrix")

#Subset counts for expressed genes and to match cellinfo 
counts <- counts[ rowSums(counts) > 0,  ]
counts <- counts[ , cellinfo$cell ]

#Normalising counts
f <- log1p(t(t(counts) / Matrix::colSums(counts))*10^6)
g <- f/log(2, base=exp(1))

#Including scale factors for zlm
cellinfo_scale <- data.frame(
  sample = cellinfo$sample,
  diagnosis = cellinfo$diagnosis,
  genes = scale(cellinfo$genes),
  age = scale(cellinfo$age),
  sex = cellinfo$sex,
  RIN = scale(cellinfo$RNA.Integrity.Number),
  PMI = scale(cellinfo$post.mortem.interval..hours.),
  region = cellinfo$region,
  Capbatch = cellinfo$Capbatch,
  Seqbatch = cellinfo$Seqbatch,
  ind = cellinfo$individual,
  ribo_perc = scale(cellinfo$RNA.ribosomal.percent))

#Create SCA object for all genes of cells in L2/3 cluster
sce_scale <- SingleCellExperiment( assays = 
                                         list( counts = g[, cellinfo$cell[ cellinfo$cluster == "L2/3"] ] ), 
                                       colData= cellinfo_scale[ cellinfo$cluster == "L2/3", ] )
sca_scale <- as(sce_scale, 'SingleCellAssay')

#Subset SCA for one gene
create_sca_gene <- function(g) {
  sca_scale <- sca_scale[ g, ]
  assay( sca_scale ) <- as.matrix( assay( sca_scale ) )
  return(sca_scale)
}

#Do zlm fitting for one gene
create_zlm_gene <- function(g) {
  sca_scale <- sca_scale[ g, ]
  assay( sca_scale ) <- as.matrix( assay( sca_scale ) )
  return(sca_scale)
  zlm_scale <- MAST::zlm( ~diagnosis + (1|ind) + genes + age + sex + RIN + PMI + 
                        region + Capbatch + Seqbatch + ribo_perc, sca_scale, method = "glmer", 
                        ebayes = F, silent=T, fitArgsD = list(nAGQ = 0))
  
  lrTest( zlm_scale, Hypothesis( "diagnosisControl" ) )
}

u <- create_sca_gene("WASH7P")

scaled_scas <- sapply( rownames(sca_scale)[1:100], create_sca_gene, simplify=FALSE )

#Permutation test on scaled
pvalue_scale <- list()

do_permutations <- function(scaclass) {
  for (i in seq(1, 100, 1)) {
    scaperm <- scaclass
    colData(scaperm) %>% as_tibble %>% left_join( by = "sample", 
       colData(scaperm) %>% as_tibble %>% dplyr::select( sample, diagnosis_old=diagnosis ) %>% 
      distinct %>% mutate( diagnosisPerm = sample(diagnosis_old) ) ) %>%
      pull( diagnosisPerm ) -> a 
    b <- colData(scaperm)
    b$diagnosisPerm <- a
    
    scaperm <- SingleCellExperiment(
      assays = list( counts =  assay(scaclass) ), colData = b)
    scaperm <- as(scaperm, 'SingleCellAssay')
    
    zlm2 <- MAST::zlm(~ diagnosisPerm + (1|ind) + genes + age + sex + RIN + PMI + 
                        region + Capbatch + Seqbatch + ribo_perc, scaperm, method = "glmer", 
                      ebayes = F, silent=T, fitArgsD = list(nAGQ = 0))
    k <- lrTest( zlm2, Hypothesis("diagnosisPermControl") )
    k <- as.data.frame(k)
    pvalue_scale[[i]] <- k$`hurdle.Pr(>Chisq)`
  }
  return(pvalue_scale)
}

p <- do_permutations(u)
  
#Looking at pvalues in QQplot
plot(-log10(ppoints(100, 0.1)), -log10(sort(unlist(p))))
abline(0,1)
