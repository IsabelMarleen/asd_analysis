#Workflow for PCA, UMAP, colsums and smoothed SATB2

#Safepoint
save.image(file="Workspace for ASD Workflow.RData")

#Setup
library( ggplot2 )
library( Matrix )
library( locfit )
library( purrr )
library( furrr )
library( tidyverse )
library(MAST)
source( "spmvar.R" )

#Load data
cellinfo <- read.delim( "~/Desktop/ASD/rawMatrix/meta.txt", stringsAsFactors=FALSE )
counts <- Seurat::Read10X("~/Desktop/ASD/rawMatrix")

#Subset counts matrix to only those cells that appear in the cellinfo table 
counts <- counts[ , cellinfo$cell ]
counts <- counts[, Matrix::colSums(counts) > 500 ]
# Remove genes that are not seen in any cell
counts <- counts[ Matrix::rowSums(counts) > 0,  ]

#Make backup of counts matrix
counts_bak2 <- counts

# Make list of sample names
samplenames <- unique( cellinfo$sample ) 

# sample info
cellinfo %>% select( sample : RNA.Integrity.Number ) %>% distinct() -> sampleTable

# Function to calculate PCA, UMAP, colSums for a sample
calc_umap_etc <- function( samplename ) {
  print( samplename )
  
  #First filter counts according to sample name
  cnts <- counts[ , cellinfo$cell[ cellinfo$sample==samplename ] ]
  
  #Select informative genes and do PCA
  frac <- t( t( cnts ) / colSums( cnts ) )
  gene_means <- rowMeans( frac )
  gene_vars <- rowVars_spm( frac )
  poisson_vmr <- mean( 1 / colSums( cnts ) )
  
  informative_genes <-  names( which( gene_vars / gene_means  >  1.5 * poisson_vmr ) )
  pca <- irlba::prcomp_irlba( t( log1p( frac[ informative_genes, ]/poisson_vmr ) ), n = 20 )$x 
  
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
add_gene("TTF2")

#Convert data into tibble for easier handling
data2 <- data %>%
  bind_rows( .id="sample" ) %>%
  as_tibble 

#save( data, file="data.rda" )

#UMAPs, with the red colour channel being the SATB2 expression
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


#UMAPs, with the red colour channel being the SATB2 expression
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


data %>%
  bind_rows( .id="sample" ) %>% 
  left_join( sampleinfo ) %>% 
  ggplot +
  geom_point( aes( UMAP1, UMAP2 , col= region), size=.3 ) +
  coord_fixed() 

#Create a plot for each sample of smoothed SATB2 against NFU1
data %>%
  bind_rows( .id="sample" ) %>%
  ggplot +
  geom_point( aes( smooth_SATB2, smooth_NFU1 ), size=.3 ) +
  facet_wrap( ~ sample )


#Using the Multimode package to predict cutoff points
#to filter for SATB2 positive cells, including all cells
#beyond the second peak line and 90% of cells between valley and second peak line

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


#Multiple testing correction for L2/3 cells
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



#Recreating MAST analysis from Schirmer paper
#https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html

# cngeneson is gene detection rate (factor recommended in MAST tutorial), 
# Capbatch is 10X capture batch, Seqbatch is sequencing batch, ind is individual label, 
# RIN is RNA integrity number, PMI is post-mortem interval and ribo_perc is 
# ribosomal RNA fraction.

MAST::zlm(~diagnosis + (1|ind) + cngeneson + age + sex + RIN + PMI + 
            region + Capbatch + Seqbatch + ribo_perc, counts, method = "glmer", 
          ebayes = F, silent=T)

f <- log1p(t(t(counts) / Matrix::colSums(counts))*10^6)
g <- f/log(2, base=exp(1))

sce <- SingleCellExperiment(assays = list(counts = g), colData= cellinfo)
sca <- as(sce, 'SingleCellAssay')
zlm <- MAST::zlm(~diagnosis + (1|individual) + genes + age + sex + RNA.Integrity.Number + 
                   post.mortem.interval..hours. + region + Capbatch + Seqbatch + 
                   RNA.ribosomal.percent, sca, method = "glmer", 
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

#Including scale factors for zlm on L2/3 and TTF2
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

sceTTF2_scale <- SingleCellExperiment( assays = 
                        list( counts = g[, cellinfo$cell[ cellinfo$cluster == "L2/3"] ] ), 
                        colData= cellinfo_scale[ cellinfo$cluster == "L2/3", ] )
scaTTF2_scale <- as(sceTTF2_scale, 'SingleCellAssay')

scaTTF2_scale <- scaTTF2_scale[ "TTF2", ]
assay( scaTTF2_scale ) <- as.matrix( assay( scaTTF2_scale ) )

zlm_scale <- MAST::zlm( ~diagnosis + (1|ind) + genes + age + sex + RIN + PMI + 
                     region + Capbatch + Seqbatch + ribo_perc, scaTTF2_scale, method = "glmer", 
                   ebayes = F, silent=T, fitArgsD = list(nAGQ = 0))

lrTest( zlm_scale, Hypothesis( "diagnosisControl" ) )


#Permutations on scaled
pvalue_scale <- list()

for (i in seq(1, 100, 1)) {
  scaTTF2perm <- scaTTF2_scale
  colData(scaTTF2perm) %>% as_tibble %>% left_join( by = "sample", 
  colData(scaTTF2perm) %>% as_tibble %>% dplyr::select( sample, diagnosis_old=diagnosis ) %>% 
                      distinct %>% mutate( diagnosisPerm = sample(diagnosis_old) ) ) %>%
    pull( diagnosisPerm ) -> a 
  b <- colData(scaTTF2perm)
  b$diagnosisPerm <- a
  
  scaTTF2perm <- SingleCellExperiment(
    assays = list( counts =  assay(scaTTF2_scale) ), colData = b)
  scaTTF2perm <- as(scaTTF2perm, 'SingleCellAssay')
  
  zlm2 <- MAST::zlm(~ diagnosisPerm + (1|ind) + genes + age + sex + RIN + PMI + 
                      region + Capbatch + Seqbatch + ribo_perc, scaTTF2perm, method = "glmer", 
                    ebayes = F, silent=T, fitArgsD = list(nAGQ = 0))
  k <- lrTest( zlm2, Hypothesis("diagnosisPermControl") )
  k <- as.data.frame(k)
  pvalue_scale[[i]] <- k$`hurdle.Pr(>Chisq)`
}

#Looking at pvalues in QQplot
plot(-log10(ppoints(100, 0.1)), -log10(sort(unlist(pvalue_scale))))

#Normal linear model
meansL23

sampleinfo2 <- cellinfo %>% 
  select( sample, diagnosis, individual, age, sex, RNA.Integrity.Number, 
          post.mortem.interval..hours., region, Capbatch, Seqbatch ) %>% 
  distinct()

fit <- limma::eBayes(limma::lmFit( meansL23[ , sampleinfo2$sample ],
              model.matrix( ~ diagnosis + age + sex + RNA.Integrity.Number + 
                              post.mortem.interval..hours. + region + Capbatch + Seqbatch , sampleinfo2 ) ))

#Permutations

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

#Permutations for DESeq
pvalue <- list()
pseudobulk_L23 <-
  sapply( sampleTable$sample, function(s)
    rowSums( counts[ , cellinfo$sample == s & cellinfo$cluster=="L2/3", drop=FALSE ] ) )

for (i in seq(1, 100, 1)) {
  columndata <- data.frame(diagnosis_new = sample(sampleTable$diagnosis))
  dds <- DESeqDataSetFromMatrix( pseudobulk_L23, columndata, ~ diagnosis_new )
  dds <- DESeq( dds )
  res <- results( dds, contrast=c("diagnosis_new", "ASD", "Control") )
  res <- res["TTF2", ]
  pvalue[[i]] <- res$pvalue
}
#padj is a list of permutations when saving the adjusted p value

plot(-log10(ppoints(100)), -log10(sort(unlist(pvalue))), 
     xlab="-log10 of evenly-spaced points", ylab="-log10 of DESeq P-Values" , main="QQ Plot DESeq")
abline(0, 1)

plot(-log10(ppoints(100)), -log10(sort(unlist(padj))), 
     xlab="-log10 of uniform Points", ylab="-log10 of P-Values" , main="Adjusted P Values")
abline(0, 1)

#Do DESeq for all different clusters and compare
pseudobulk_L23 <-
  sapply( sampleTable$sample, function(s)
    rowSums( counts[ , cellinfo$sample == s & cellinfo$cluster=="L2/3", drop=FALSE ] ) )
dds_L23 <- DESeqDataSetFromMatrix( pseudobulk_L23, sampleTable, ~ diagnosis )
dds_L23 <- DESeq( dds_L23 )
res_L23 <- results( dds_L23 )

sapply( unique( cellinfo$cluster ), function( clust ) {
    assign( paste0( "pseudobulk_", clust ),  sapply( sampleTable$sample, function(s) {
           rowSums( counts[ , cellinfo$sample == s & cellinfo$cluster == clust, drop=FALSE ] ) }) )
  
    pseudobulk <-  sapply( sampleTable$sample, function(s) {
    rowSums( counts[ , cellinfo$sample == s & cellinfo$cluster == clust, drop=FALSE ] ) } )
    dds <- DESeqDataSetFromMatrix( pseudobulk, sampleTable, ~ diagnosis )
    dds <-  DESeq( dds )
    assign( paste0( "res_", clust ), results( dds) )

}  )

cluster <- "Microglia"
do_DESeq_cluster <- function( cluster ){
  pseudobulk <- sapply( sampleTable$sample, function(s)
      rowSums( counts[ , cellinfo$sample == s & cellinfo$cluster==cluster, drop=FALSE ] ) )
  dds <- DESeqDataSetFromMatrix( pseudobulk, sampleTable, ~ diagnosis)
  keep <- rowSums(counts(dds)) >= 10
  keep2 <- colSums(counts(dds)) > 0
  dds <- dds[keep,keep2]
  dds <- DESeq( dds )
  return(dds)
}


dds_L23 <- do_DESeq_cluster("L2/3")
dds_L4 <- do_DESeq_cluster("L4")
dds_L56 <- do_DESeq_cluster("L5/6")
dds_L56CC <- do_DESeq_cluster("L5/6-CC")
dds_NRGNII <- do_DESeq_cluster("Neu-NRGN-II")
dds_NRGNI <- do_DESeq_cluster("Neu-NRGN-I")
dds_MG <- do_DESeq_cluster("Microglia")
dds_OD <- do_DESeq_cluster("Oligodendrocytes")
dds_OPC <- do_DESeq_cluster("OPC")
dds_ASTFB <- do_DESeq_cluster("AST-FB")
dds_ASTPP <- do_DESeq_cluster("AST-PP")
dds_Nmat <- do_DESeq_cluster("Neu-mat")
dds_INSST <- do_DESeq_cluster("IN-SST")
dds_INPV <- do_DESeq_cluster("IN-PV")
dds_INVIP <- do_DESeq_cluster("IN-VIP")
dds_INSV2C <- do_DESeq_cluster("IN-SV2C")
dds_ET <- do_DESeq_cluster("Endothelial")


#Comparing to genes from paper
results_lucas <- read.delim(file="~/Desktop/ASD/S4.csv", sep=";", dec=",")

plot_cluster <- function(clust1, dds) {
  left_join(as_tibble(filter(results_lucas, results_lucas$X...Cell.type==clust1)), 
            as_tibble(results(dds), rownames="Gene.name")) %>%
    ggplot+
    geom_point(aes(log(q.value), log(padj)))+
    geom_vline(aes(xintercept=-1))+
    geom_hline(aes(yintercept=-1))+
    ggtitle(clust1)
}

plot_cluster("L2/3", dds_L23)
plot_cluster("L4", dds_L4)
plot_cluster("L5/6", dds_L56)
plot_cluster("L5/6-CC", dds_L56CC)
plot_cluster("AST-FB", dds_ASTFB)
plot_cluster("AST-PP", dds_ASTPP)
plot_cluster("Endothelial", dds_ET)
plot_cluster("IN-PV", dds_INPV)
plot_cluster("IN-SST", dds_INSST)
plot_cluster("IN-SV2C", dds_INSV2C)
plot_cluster("IN-VIP", dds_INVIP)
plot_cluster("Microglia", dds_MG)
plot_cluster("Neu-NRGN-I", dds_NRGNI)
plot_cluster("Neu-NRGN-II", dds_NRGNII)
plot_cluster("Neu-mat", dds_Nmat)
plot_cluster("Oligodendrocytes", dds_OD)
plot_cluster("OPC", dds_OPC)


#Look at differential gene expression
plot_hist_gene <- function(gene) {
  ans <- left_join( data2, select( sampleTable, diagnosis, sample) )%>%
    dplyr::filter(cellinfo$cluster == "L2/3")%>%
    mutate(state=case_when(
      diagnosis == "ASD" ~ "asd",
      TRUE ~ "control" )) %>%
    dplyr::filter( ., ans[[ paste0( "smooth_", gene ) ]] < .9 )
  ggplot(ans)+
    geom_density(aes(ans[[paste0("smooth_", gene)]]^.15, col=state, group=ans$sample))
}
plot_hist_gene_celltypes("TTF2")

plot_hist_gene2 <- function(gene) {
  ans <- left_join( data2, select( sampleTable, diagnosis, sample) )%>%
    dplyr::filter(cellinfo$cluster == "L2/3")%>%
    mutate(state=case_when(
      diagnosis == "ASD" ~ "asd",
      TRUE ~ "control" )) %>%
    dplyr::filter( ., ans[[ paste0( "smooth_", gene ) ]] < .9 )
  ggplot(ans)+
    geom_density(aes(ans[[paste0("smooth_", gene)]]^.15, col=state))+
    facet_wrap(~sample)
}
plot_hist_gene2("TTF2")




##Create new clusters by using louvain--------------------------------------------------
library(igraph)
a <- data2 %>%
    select(2:21) %>%
    as.matrix()%>%
    FNN::get.knn()
g <- graph_from_edgelist( as.matrix( map_dfr( 1:10, function(i)
  data.frame( from=1:nrow(a$nn.index), to=a$nn.index[,i] ) ) ) )
louv <- cluster_louvain( as.undirected( g, "collapse" ) )
nn <- membership(louv)

data2%>%
  filter(sample=="1823_BA24")%>%
  {
    ggplot(.)+
      geom_point(aes(.$UMAP1, .$UMAP2, col=nn))
  }

