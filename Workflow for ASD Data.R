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
cellinfo <- read.delim( "~/Desktop/rawMatrix/meta.txt", stringsAsFactors=FALSE )
counts <- Seurat::Read10X("~/Desktop/rawMatrix")

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

save( data, file="data.rda" )

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
       green = pmin( 1, a$smooth_COX5B^.15/.3),
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


#Using the MAST package in an attempt to replicate Schirmer results
MAST::zlm(~diagnosis + (1|ind) + cngeneson + age + sex + RIN + PMI + 
      region + Capbatch + Seqbatch + ribo_perc, counts, method = "glmer", 
    ebayes = F, silent=T)

# Where cngeneson is gene detection rate (factor recommended in MAST tutorial), 
# Capbatch is 10X capture batch, Seqbatch is sequencing batch, ind is individual label, 
# RIN is RNA integrity number, PMI is post-mortem interval and ribo_perc is 
# ribosomal RNA fraction.

#Following protocol
#https://www.bioconductor.org/packages/release/bioc/vignettes/MAST
#/inst/doc/MAITAnalysis.html


counts3 <- counts

freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)
data(maits, package="MAST")
dim(maits$expressionmat)
head(maits$cdat)
head(maits$fdat)

scaRaw <- MAST::FromMatrix(t(maits$expressionmat), maits$cdat, maits$fdat)

f <- log1p(t(t(counts) / colSums(counts))*10^6) #%>%
g <- f/log(2, base=exp(1))
