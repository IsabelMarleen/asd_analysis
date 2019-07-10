#Workflow for PCA, UMAP, colsums and smoothed SATB2

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

#Subset counts matrix to onlz those cells that appear in the cellinfo table 
counts <- counts[ , cellinfo$cell ]

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

#Create a plot for each sample of smoothed SATB2 against NFU1
data %>%
  bind_rows( .id="sample" ) %>%
  ggplot +
  geom_point( aes( smooth_SATB2, smooth_NFU1 ), size=.3 ) +
  facet_wrap( ~ sample )
