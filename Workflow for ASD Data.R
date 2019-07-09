#Workflow for PCA, UMAP, colsums and smoothed SATB2

#Setup
library( ggplot2 )
library( Matrix )
library( locfit )
library( purrr )
library( furrr )
library( tidyverse )
source( "spmvar.R" )
cellinfo <- read.delim("~/Desktop/rawMatrix/meta.txt")
counts <- Seurat::Read10X("~/Desktop/rawMatrix")
counts_bak <- counts

##Make list of sample names
samplenames <- as.data.frame(select(cellinfo, sample), stringsasFactors=FALSE, drop=FALSE )  %>%
  distinct()#%>%
  #pull(sample)

#Filtering counts
counts <- counts[ , colSums( counts>0 ) >= 1000 ]
counts <- counts[ rowSums(counts) > 0,  ]

#Function
analyse <- function (samplename, counts) {
  
  ans = list()

  #First filter counts according to sample name, live cells and exressed genes
  counts <- counts[,as.character(cellinfo$cell[cellinfo$sample==samplename])]
  
  #Select informative genes and do PCA
  frac <- t(t(counts) / colSums(counts))
  gene_means <- rowMeans( frac )
  gene_vars <- rowVars_spm( frac )
  poisson_vmr <- mean( 1 / colSums( counts ) )
  
  informative_genes <-  names(which(gene_vars / gene_means  >  1.5 * poisson_vmr ))
  ans$pca <- irlba::prcomp_irlba( t(log1p(frac[informative_genes,]/poisson_vmr)), n = 20)$x
  
  #Make colsums list
  ans$colsums <- colSums(counts)
  
  #Do UMAP
  
  ans$umap <- uwot::umap( asd$pca,n_neighbors = 30, min_dist = .3, metric = "cosine" )
  
  #Smooth SATB2
  ans$rawSATB2 <- counts["SATB2",]
  ans$smoothedSATB2 <-  suppressWarnings(predict( locfit.raw( 
    ans$pca[,1:15], counts["SATB2",], base = log( colSums(counts) ),
    alpha = .1, deg = 1, family = "poisson", ev = dat() ) ) ) 
  
  #Show smoothed UMAP
  # return(ggplot( NULL, aes( x=asd$umap[,1], y=asd$umap[,2] ) ) + coord_fixed() +
  #   geom_point( aes( col = asd$smoothedSATB2^.15 ), size=.15 ) + labs( col="SATB2") + 
  #   scale_color_gradientn( colours = rje::cubeHelix(100)[1:80], limits = c( 0.12, .45 ), 
  #                          breaks = (tics*xi)^.15, labels=tics, oob=scales::squish))  
  # 
  # 
  #Align with Clustering from Paper
 # return(data.frame(cell=colnames(counts), asd$umap) %>%
  #  left_join(cellinfo, by = "cell") %>%
  #  ggplot + geom_point(aes(x=X1, y=X2, colour=cluster), size=0.15) + coord_fixed())
  #return(plotly::ggplotly())
  return(ans)
}

#Command that works through first two samples
x <- sapply( samplenames$sample[2], analyse, counts, simplify=FALSE )

#Command that assigns samplenames to sublists
names(x) <- samplenames$sample

#back-up of x
x_bak <- x

#Command that works through all samples
x <- sapply( samplenames$sample, analyse, counts, simplify=FALSE )

#Makes tibble out of list, then feeds into ggplot to produce all umaps
xz <- map_dfr( x, .id="sample", function(xx) 
  tibble( colSums=xx$colsums, umapX=xx$umap[,1], umapY=xx$umap[,2], 
     rawSATB2=xx$rawSATB2, smoothedSATB2=xx$smoothedSATB2) ) 
xz %>%
ggplot +
  geom_point( aes( x=umapX, y=umapY, col=smoothedSATB2^.15 ), size=.1 ) +
  coord_fixed() +
  facet_wrap( ~ sample )+
  scale_color_gradientn( colours = rje::cubeHelix(100)[1:80], limits = c( 0.12, .45 ), 
          breaks = (tics*xi)^.15, labels=tics, oob=scales::squish)


##UMAPS for NFU1
samplename <- "1823_BA24"

analyse_NFU1 <- function (samplename, counts) {
  
  ans = list()
  
  #First filter counts according to sample name, live cells and exressed genes
  #Think about using dplyr::filter and piping instead
  counts <- counts[, as.character(cellinfo$cell[cellinfo$sample==samplename])]
  
  #Select informative genes and do PCA
  frac <- t(t(counts) / colSums(counts))
  gene_means <- rowMeans( frac )
  gene_vars <- rowVars_spm( frac )
  poisson_vmr <- mean( 1 / colSums( counts ) )
  
  informative_genes <-  names(which(gene_vars / gene_means  >  1.5 * poisson_vmr ))
  ans$pca <- irlba::prcomp_irlba( t(log1p(frac[informative_genes,]/poisson_vmr)), n = 20)$x
  
  #Make colsums list
  ans$colsums <- colSums(counts)
  
  #Do UMAP
  
  ans$umap <- uwot::umap( ans$pca,n_neighbors = 30, min_dist = .3, metric = "cosine" )
  
  #Smooth SATB2
  ans$rawSATB2 <- counts["SATB2",]
  ans$smoothedSATB2 <-  suppressWarnings(predict( locfit.raw( 
    ans$pca[,1:15], counts["SATB2",], base = log( colSums(counts) ),
    alpha = .1, deg = 1, family = "poisson", ev = dat() ) ) ) 
  
  #Smooth NFU1
  ans$rawNFU1 <- counts["NFU1",]
  ans$smoothedNFU1 <-  suppressWarnings(predict( locfit.raw( 
    ans$pca[,1:15], counts["NFU1",], base = log( colSums(counts) ),
    alpha = .1, deg = 1, family = "poisson", ev = dat() ) ) ) 
  
  #Show smoothed UMAP
  # return(ggplot( NULL, aes( x=asd$umap[,1], y=asd$umap[,2] ) ) + coord_fixed() +
  #   geom_point( aes( col = asd$smoothedSATB2^.15 ), size=.15 ) + labs( col="SATB2") + 
  #   scale_color_gradientn( colours = rje::cubeHelix(100)[1:80], limits = c( 0.12, .45 ), 
  #                          breaks = (tics*xi)^.15, labels=tics, oob=scales::squish))  
  # 
  # 
  #Align with Clustering from Paper
  # return(data.frame(cell=colnames(counts), asd$umap) %>%
  #  left_join(cellinfo, by = "cell") %>%
  #  ggplot + geom_point(aes(x=X1, y=X2, colour=cluster), size=0.15) + coord_fixed())
  #return(plotly::ggplotly())
  return(ans)
}

#First collect NFU1 data
k <- sapply( samplenames$sample, analyse_NFU1, counts, simplify=FALSE )

xk$rawNFU1 <- counts["NFU1",]
#xx$smoothedNFU1 <-  suppressWarnings(predict( locfit.raw( 
 # x$pca[,1:15], counts["NFU1",], base = log( colSums(counts) ),
  #alpha = .1, deg = 1, family = "poisson", ev = dat() ) ) )

map_dfr( k, .id="sample", function(xx) 
  tibble(umapX=xx$umap[,1], umapY=xx$umap[,2], 
          smoothedNFU1=xx$smoothedNFU1) ) %>%
  ggplot +
  geom_point( aes( x=umapX, y=umapY, col=smoothedNFU1^.15 ), size=.1 ) +
  coord_fixed() +
  facet_wrap( ~ sample )+
  scale_color_gradientn( colours = rje::cubeHelix(100)[1:80], limits = c( 0.12, .45 ), 
                         breaks = (tics*xi)^.15, labels=tics, oob=scales::squish)


  
#Using the Multimode package to predict where the lines need to go
a <- sapply(names(x), function (s) {
  a <- x[[s]]$smoothedSATB2
  a <- a[a<.9]
  multimode::locmodes(a^.15, 2 )$location
})

#b is a tibble including the positions of the three lines for each sample
b <- a %>% t %>% as_tibble( rownames="sample" ) %>%
  gather( mode, pos, V1:V3 )

#UMAPS for SATB2
xz %>%
  filter( smoothedSATB2^.15 < .9) %>%
  ggplot(aes( x=x[[s]]^.15 ) ) +
  facet_wrap( ~ sample, scales = "free" ) 


ans$umap <- uwot::umap( asd$pca,n_neighbors = 30, min_dist = .3, metric = "cosine" )


#Histograms for SATB2
xz %>%
  filter( smoothedSATB2^.15 < .9) %>%
  ggplot(aes( x=smoothedSATB2^.15 ) ) +
  geom_histogram() +
  facet_wrap( ~ sample, scales = "free" ) +
  geom_vline( data=b, aes( xintercept=pos ), color="yellow")


xz %>%
  group_by(sample) %>%
  filter(smoothedSATB2^.15 < .9) %>%
  summarise(modlocs=list(multimode::locmodes(smoothedSATB2^.15,2)))


#Histogram of "5163_BA24" for SATB2
xz%>%
  filter(sample=="5163_BA24") %>%
  filter(smoothedSATB2^.15 < .9) %>%
  ggplot(aes( x=smoothedSATB2^.15 ) ) +
  geom_histogram() +
  facet_wrap( ~ sample, scales = "free" ) +
  geom_vline( data=c, aes( xintercept=pos ), color="yellow")

c <- filter(b, sample=="5163_BA24")


#Use predicted lines by locfit to filter for SATB2 positive cells, including all cells
#beyond the second peak line and 90% of cells between valley and second peak line
#Use function filter_SATB2, which takes raw data and a cutoff point

filter_SATB2(counts, b$mode["V3"])

#filter2_SATB2 <- function(raw, filter1, filter2) {
#   SATB2pos <- as.data.frame(t(raw)) %>%
#     rownames_to_column("cell") %>%
#     filter(SATB2>=filter1) %>% #filter1 = V3
#     filter(SATB2>=(filter2-filter1)*.1) #filter2=V2
#     
#     pull(cell)
#   
#   return(SATB2pos)
# }




# <- filter2_SATB2(counts, b$mode["V3"], b$mode["V2"])

# kk <- map_dfr( x, .id="sample", function(xx) 
#   tibble(cell=cellinfo$cell, smoothedSATB2=xx$smoothedSATB2) ) 
# 


for( s in names(x) ) {
  a <- x[[s]]$smoothedSATB2
  a <- a[ a < .9 ]
  x[[s]]$locmodes <- multimode::locmodes( a ^.15, 2 )$location^(1/.15)
  x[[s]]$SATB2thresh <- x[[s]]$locmodes[2]
}

#Making sure that cellinfo has the same cells as counts
cellinfo <- cellinfo[as.character(cellinfo$cell) %in% colnames(counts),]

s <- "1823_BA24"
g <- "RPL18"

avg_genes_for_SATB2pos <- function( s ) {
  #Subsetting names of raw SATB2 with smoothed SATB2 above threshhold
  SATB2pos <- names(x[[s]]$rawSATB2)[ x[[s]]$smoothedSATB2 > x[[s]]$SATB2thresh ]
  # Get fractions for these, for gene g, and take average
  rowMeans( t( t(counts[ , SATB2pos ]) / colSums( counts[ , SATB2pos ] ) ) )
}

means_SATB2pos <- sapply( names(x), avg_genes_for_SATB2pos )

diagnoses <- cellinfo %>% select( sample, diagnosis ) %>% distinct() %>% deframe()

stopifnot( all( names(diagnoses) == colnames(means_SATB2pos) ) )

ttres <- genefilter::rowttests( means_SATB2pos, diagnoses )

ggplot( ttres ) +
  geom_point( aes( x=dm, y=-log10(p.value) ) ) +
  scale_x_continuous( limits = c( -0.002, 0.002 ) )

tibble( NFU1=means_SATB2pos["NFU1",], diagnoses ) %>% ggplot + geom_point(aes(x=1,y=NFU1, diagnoses=))
