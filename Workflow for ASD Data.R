#Workflow for PCA, UMAP, colsums and smoothed SATB2

#Safepoint
#save.image(file="Workspace for ASD Workflow.RData")

#Setup
library( ggplot2 )
library( Matrix )
library( locfit )
library( purrr )
library( furrr )
library( tidyverse )
source( "spmvar.R" )

#Load data
cellinfo <- read.delim( "~/sds/sd17l002/u/isabel/rawMatrix/meta.txt", stringsAsFactors=FALSE )
counts <- Seurat::Read10X("~/sds/sd17l002/u/isabel/rawMatrix")

#Subset counts matrix to only those cells that appear in the cellinfo table 
counts <- counts[ , cellinfo$cell ]
counts <- counts[, colSums(counts) > 500 ]
# Remove genes that are not seen in any cell
counts <- counts[ rowSums(counts) > 0,  ]

#Make backup of counts matrix
counts_bak2 <- counts

# Make list of sample names
samplenames <- unique( cellinfo$sample ) 

# sample info
cellinfo %>% select( sample : RNA.Integrity.Number ) %>% distinct() -> sampleinfo

# Function to calculate PCA, UMAP, colSums for a sample
calc_umap_etc <- function( samplename ) {
  print( samplename )
  
  #First filter counts according to sample name
  cnts <- counts[ , cellinfo$cell[ cellinfo$sample==samplename ] ]
  
  #Select informative genes and do PCA
  frac <- t(t(cnts) / colSums(cnts))
  gene_means <- rowMeans( frac )
  gene_vars <- rowVars_spm( frac )
  poisson_vmr <- mean( 1 / colSums( cnts ) )
  
  informative_genes <-  names(which(gene_vars / gene_means  >  1.5 * poisson_vmr ))
  pca <- irlba::prcomp_irlba( t(log1p(frac[informative_genes,]/poisson_vmr)), n = 20)$x 
  
  #Do UMAP
  umap <- uwot::umap( pca, n_neighbors = 30, min_dist = .3, metric = "cosine" )
  colnames(umap) <- c( "UMAP1", "UMAP2" )
  
  # Make data frame, add two more useful columns
  ans <- cbind( as.data.frame(pca), as.data.frame(umap) )
  
  ans$cluster <- cellinfo$cluster[ cellinfo$sample==samplename ]
  ans$colsums <- colSums(cnts)
  rownames(ans) <- cellinfo$cell[ cellinfo$sample==samplename ]
  
  ans  
}  

#Calculate UMAP, PCA etc for all samples  
data <- sapply( samplenames, calc_umap_etc, simplify=FALSE )



#Function that adds raw and smoothed gene expression to data
add_gene <- function( gene )
{
  for( samplename in names(data) ) {
    cat( samplename, " " )
    data[[ samplename ]][[ paste0( "raw_", gene ) ]] <<- counts[ gene, rownames(data[[samplename]]) ]
    data[[ samplename ]][[ paste0( "smooth_", gene ) ]] <<- 
       suppressWarnings(predict( locfit.raw( 
         as.matrix( data[[ samplename ]][ ,paste0( "PC", 1:15 ) ] ), 
         data[[ samplename ]][[ paste0( "raw_", gene ) ]],
         base = log( data[[ samplename ]]$colsums ),
         alpha = .1, deg = 1, family = "poisson", ev = dat() ) ) ) 
  }
}
  
add_gene( "SATB2" )
add_gene( "NFU1" )
add_gene("COX5B")

#save( data, file="data.rda" )

#Create all UMAPs, whith the red colour channel being the SATB2 expression
#and the green one the COXB5 expression, overlapping regions are yellow
#No systematic difference between ASD and Control samples visible
data %>%
bind_rows( .id="sample" ) %>% 
left_join( sampleinfo ) %>%
unite( smpl, individual, region, diagnosis ) -> a
a %>% ggplot +
  geom_point( aes( UMAP1, UMAP2 ), size=.3,
     col = rgb( 
       red = pmin( 1, a$smooth_SATB2^.15/.4), 
       green = pmin( 1, a$smooth_NFU1^.15/.3),
       blue = .3 ) ) +
  coord_fixed() +
  facet_wrap( ~ smpl ) +
  theme_dark() + theme( plot.background = element_rect(fill="black") )


#Create all UMAPs, whith the red colour channel being the SATB2 expression
#and the green one the NFU1 expression, overlapping regions are yellow
data %>%
  bind_rows( .id="sample" ) %>% 
  left_join( sampleinfo ) %>%
  unite( smpl, individual, region, diagnosis ) -> a
a %>% ggplot +
  geom_point( aes( UMAP1, UMAP2 ), size=.3,
              col = rgb( 
                red = pmin( 1, a$smooth_SATB2^.15/.4), 
                green = pmin( 1, a$smooth_NFU1^.15/.3),
                blue = .3 ) ) +
  coord_fixed() +
  facet_wrap( ~ smpl ) +
  theme_dark() + theme( plot.background = element_rect(fill="black") )

#Create a plot for each sample of smoothed SATB2 against NFU1
data %>%
  bind_rows( .id="sample" ) %>%
  ggplot +
  geom_point( aes( smooth_SATB2, smooth_NFU1 ), size=.3 ) +
  facet_wrap( ~ sample )



#Not cleaned yet



#Using the Multimode package to predict where the lines in the histograms need to go
lines <- sapply( names(data), function (s) {
  lines <- data[[s]]$smooth_SATB2
  lines <- lines[lines<.9]
  multimode::locmodes(lines^.15, 2 )$location
})

#modes_positions is a tibble including the positions of the three lines for each sample
modes_positions <- lines %>% t %>% as_tibble( rownames="sample" ) %>%
  gather( mode, pos, V1:V3 )



#Use predicted lines by locfit to filter for SATB2 positive cells, including all cells
#beyond the second peak line and 90% of cells between valley and second peak line

#Creating new list
location <- list()
for( s in names(data) ) {
  a <- data[[s]]$smooth_SATB2
  a <- a[ a < .9 ]
  location[[s]]$locmodes <- multimode::locmodes( a ^.15, 2 )$location^(1/.15)
  location[[s]]$SATB2thresh <- location[[s]]$locmodes[2]
}


##Averaging genes in SATB2pos cells to compare ASD/control using t tests
#Making sure that cellinfo has the same cells as counts
cellinfo <- cellinfo[(cellinfo$cell) %in% colnames(counts),]

avg_genes_for_SATB2pos <- function( s ) {
  #Subsetting names of smooth SATB2 with rownames above threshhold
  SATB2pos <- rownames(data[[s]])[ data[[s]]$smooth_SATB2 > location[[s]]$SATB2thresh ]
  # Get fractions for these, for gene g, and take average
  rowMeans( t( t(counts[ , SATB2pos ]) / colSums( counts[ , SATB2pos ] ) ) )
}

means_SATB2pos <- sapply( names(data), avg_genes_for_SATB2pos )

diagnoses <- cellinfo %>% select( sample, diagnosis ) %>% distinct() %>% deframe() %>%
  as.factor()

stopifnot( all( names(diagnoses) == colnames(means_SATB2pos) ) )

#t Test
ttres <- genefilter::rowttests( means_SATB2pos, diagnoses )
rownames(ttres) <- rownames(means_SATB2pos)

#Multiple testing correction
ttres$padj <- p.adjust(ttres$p.value, method="BH")

ggplot( ttres ) +
  geom_point( aes( x=log10(dm), y=-log10(p.value) ) ) +
  scale_x_continuous( limits = c( -0.002, 0.002 ) )

tibble( NFU1=means_SATB2pos["NFU1",], diagnoses ) %>% ggplot + geom_point(aes(x=1,y=NFU1, diagnoses=))

#Replacing zero-inflated model with simple t-tests on all cells, not only SATB2 positive ones

avg_genes <- function( s ) {
  #Subsetting names of data according to sample
  allcells <- rownames(data[[s]])
  # Get fractions for gene g, and take average
  rowMeans( t( t(counts[ , allcells ]) / colSums( counts[ , allcells] ) ) )
}

means <- sapply( names(data), avg_genes)

#Assertion
stopifnot( all( names(diagnoses) == colnames(means) ) )

#t Test for all cells
ttres2 <- genefilter::rowttests( means, diagnoses )
rownames(ttres2) <- rownames(means)


#Multiple testing correction for all cells
#Still no siginificant evidence for differences
ttres2$padj <- p.adjust(ttres2$p.value, method="BH")


#Using the MAST package in an attempt to replicate Schirmer results
MAST::zlm(~diagnosis + (1|ind) + cngeneson + age + sex + RIN + PMI + 
      region + Capbatch + Seqbatch + ribo_perc, counts, method = "glmer", 
    ebayes = F, silent=T)

#t Test for only L2/3
avg_genesL23 <- function( s ) {
  SATB2pos <- rownames(data[[s]])
  #Subsetting names of data according to sample and cluster
  allcells <- rownames(data[[s]])[ data[[s]]$cluster == "L2/3" ]
  # Get fractions for gene g, and take average
  rowMeans( t( t(counts[ , allcells ]) / colSums( counts[ , allcells] ) ) )
}

meansL23 <- sapply( names(data), avg_genes)

#Assertion
stopifnot( all( names(diagnoses) == colnames(meansL23) ) )

#t Test for all L2/3 cells
ttres3 <- genefilter::rowttests( meansL23, diagnoses )
rownames(ttres3) <- rownames(meansL23)


#Multiple testing correction for all cells
#Still no siginificant evidence for differences
ttres3$padj <- p.adjust(ttres3$p.value, method="BH")

#Relating results of t test with results given in paper for LMM
S4 <- read.table("~/Desktop/S4.csv", header= TRUE, sep=";")
S4 <- S4[S4$X...Cell.type == "L2/3", ]
S4$Gene.name %in% rownames(ttres3)


#Plot difference in means against p Value and colour in genes from S4
ggplot() +
  geom_point(aes(x=ttres3$dm, y=-log10(ttres3$p.value), col=rownames(ttres3)%in% S4$Gene.name)) +
  scale_x_continuous(limits=c(-0.001,0.001),oob=scales::squish)
  
 


#Plot ttres against ttres2 to see whether filtering for SATB2 pos cells makes a systematic difference
plot(ttres$p.value[1:100], ttres2$p.value[1:100], type="p")
plot(ttres$padj[1:100], ttres2$padj[1:100], type="p")



#Following protocol
#https://www.bioconductor.org/packages/release/bioc/vignettes/MAST
#/inst/doc/MAITAnalysis.html

# cngeneson is gene detection rate (factor recommended in MAST tutorial), 
# Capbatch is 10X capture batch, Seqbatch is sequencing batch, ind is individual label, 
# RIN is RNA integrity number, PMI is post-mortem interval and ribo_perc is 
# ribosomal RNA fraction.

library(MAST)

counts3 <- counts

f <- log1p(t(t(counts) / Matrix::colSums(counts))*10^6)
g <- f/log(2, base=exp(1))

sce <- SingleCellExperiment(assays = list(counts = g), colData= cellinfo)
sca <- as(sce, 'SingleCellAssay')
zlm <- MAST::zlm(~diagnosis + (1|individual) + genes + age + sex + RNA.Integrity.Number + post.mortem.interval..hours. + 
                   region + Capbatch + Seqbatch + RNA.ribosomal.percent, sca, method = "glmer", 
                 ebayes = F, silent=T)

#Fitting model only for NFU1

zlm <- MAST::zlm(~diagnosis + (1|individual) + genes + age + sex + RNA.Integrity.Number + post.mortem.interval..hours. + 
                   region + Capbatch + Seqbatch + RNA.ribosomal.percent, scaNFU1, method = "glmer", 
                 ebayes = F, silent=T)

lrTest( zlm, Hypothesis("diagnosisControl") )


#Fitting model with only L2/3 cells and TTF2

sceTTF2 <- SingleCellExperiment( assays = 
                                list( counts = g[, cellinfo$cell[ cellinfo$cluster == "L2/3" ] ] ), 
                                colData= cellinfo[ cellinfo$cluster == "L2/3", ] )
scaTTF2 <- sca[ "TTF2", ]
assay( scaTTF2 ) <- as.matrix( assay( scaTTF2 ) )

zlm2 <- MAST::zlm( ~diagnosis + (1|individual) + genes + age + sex + RNA.Integrity.Number + post.mortem.interval..hours. + 
                    region + Capbatch + Seqbatch + RNA.ribosomal.percent, scaTTF2, method = "glmer", 
                  ebayes = F, silent=T, fitArgsD = list(nAGQ = 0))

lrTest( zlm2, Hypothesis( "diagnosisControl" ) )

#Normal linear model
meansL23

sampleinfo2 <- cellinfo %>% select( sample, diagnosis, individual, age, sex, RNA.Integrity.Number, 
                                    post.mortem.interval..hours., region, Capbatch, Seqbatch ) %>% distinct()

fit <- limma::eBayes(limma::lmFit( meansL23[ , sampleinfo2$sample ],
                                   model.matrix( ~ diagnosis + age + sex + RNA.Integrity.Number + post.mortem.interval..hours. + 
                                                   region + Capbatch + Seqbatch , sampleinfo2 ) ))


#Permutations test

pvalue <- list()

for (i in seq(1, 10, 1)) {
  scaTTF2perm <- scaTTF2
  colData(scaTTF2perm) %>% as_tibble %>% left_join( by = "sample", 
    colData(scaTTF2perm) %>% as_tibble %>% dplyr::select( sample, diagnosis_old=diagnosis ) %>% 
    distinct %>% mutate( diagnosisPerm = sample(diagnosis_old) ) ) %>%
    pull( diagnosisPerm ) -> a 
  b <- colData(scaTTF2perm)
  b$diagnosisPerm <- a
  
  scaTTF2perm <- SingleCellExperiment(
    assays = list( counts =  assay(scaTTF2) ), colData = b)
  scaTTF2perm <- as(scaTTF2perm, 'SingleCellAssay')
  
  zlm2 <- MAST::zlm(~ diagnosisPerm + (1|individual) + genes + age + sex + RNA.Integrity.Number + post.mortem.interval..hours. + 
                      region + Capbatch + Seqbatch + RNA.ribosomal.percent, scaTTF2perm, method = "glmer", 
                    ebayes = F, silent=T, fitArgsD = list(nAGQ = 0))
 k <- lrTest( zlm2, Hypothesis("diagnosisPermControl") )
 k <- as.data.frame(k)
  pvalue[[i]] <- k$`hurdle.Pr(>Chisq)`
}




scaTTF2perm <- scaTTF2
colData(scaTTF2perm) %>% as_tibble %>% left_join( by = "sample", 
   colData(scaTTF2perm) %>% as_tibble %>% dplyr::select( sample, diagnosis_old=diagnosis ) %>% 
      distinct %>% mutate( diagnosisPerm = sample(diagnosis_old) ) ) %>%
  pull( diagnosisPerm ) -> a 
b <- colData(scaTTF2perm)
b$diagnosisPerm <- a
colData(scaTTF2perm) <- b

scaTTF2perm <- SingleCellExperiment(
  assays = list( counts =  assay(scaTTF2) ), colData = b)
scaTTF2perm <- as(scaTTF2perm, 'SingleCellAssay')

zlm2 <- MAST::zlm(~ diagnosis + (1|individual) + genes + age + sex + RNA.Integrity.Number + post.mortem.interval..hours. + 
                    region + Capbatch + Seqbatch + RNA.ribosomal.percent, scaTTF2, method = "glmer", 
                  ebayes = F, silent=T, fitArgsD = list(nAGQ = 0))
lrTest( zlm2, Hypothesis("diagnosisControl") )








   