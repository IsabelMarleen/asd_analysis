# This script crashes with 16 GB or less of RAM.
# Rsession will use 30 GB or RAM in the long-run, not sure about peaks.





# Load --------------------------------------------------------------------


library( tidyverse )
library( Matrix )
library( irlba )
library( uwot )
library( FNN )
library( igraph )
library( cowplot )

# convenience functions such as col_pwr_trans, rowVars_spm etc.:
scr_dir <- "/home/frauhammer/sc_methods_dev/src/"
source(file.path(scr_dir, "functions_universal.R"))


path <- "/home/frauhammer/sds_copy/ASD/"
cellinfo <- read.delim( file.path( path, "meta.txt" ), stringsAsFactors=FALSE )
counts <- readMM( file.path( path, "matrix.mtx" ) )
gene_info <- read.delim( file.path( path, "genes.tsv" ), header=FALSE, as.is=TRUE ) %>%
  mutate(unique = case_when(
  duplicated(V2) | duplicated(V2, fromLast=T) ~ paste(V2, V1, sep="_"),
  TRUE ~ V2))
rownames(counts) <- gene_info$unique
colnames(counts) <- readLines( file.path( path, "barcodes.tsv" ) )


sampleTable <-
  cellinfo %>% select( sample : RNA.Integrity.Number ) %>% unique
sampleTable

# extracting gene expression is much faster in column-sparse format:
Tcounts <- as(t(counts), "dgCMatrix")   #  fast:   Tcounts[, "SYN1"]

Ccounts <- as(counts, "dgCMatrix")      #  fast:   Ccounts[, 1337]   and  colSums(Ccounts)







# Preprocessing -----------------------------------------------------------

# load (or re-execute everything in this section):
sfs <- colSums(counts)
norm_counts <- t(t(Ccounts) / colSums(Ccounts))
rownames(norm_counts) <- rownames(Ccounts)
load(file.path("~", "asd_analysis", "savepoint", "umap_euc_spread10.RData"))




# informative genes, PCA, UMAP:
poisson_vmr <- mean(1/sfs)
gene_means <- rowMeans( norm_counts )
gene_vars <- rowVars_spm( norm_counts )
cells_expressing <- rowSums( counts != 0 )
is_informative <- gene_vars/gene_means > 1.5 * poisson_vmr  &  cells_expressing > 100
plot(gene_means, gene_vars/gene_means, pch=".", log = "xy")
points(gene_means[is_informative], (gene_vars/gene_means)[is_informative], pch=".", col = "red" )
pca <- irlba::prcomp_irlba( x = sqrt(t(norm_counts[is_informative,])),
                            n = 40,
                            scale. = TRUE)
umap_euc <- uwot::umap( pca$x, spread = 10, n_threads = 40)
umap_cos <- uwot::umap( pca$x, metric = "cosine", spread = 10, n_threads = 40)

# save(umap_euc,
#      file = file.path("~", "asd_analysis", "savepoint", "umap_euc_spread10.RData"))






# Clusters ---------------------------------------------------


# load (or re-execute everything in this section):
load(file.path("~", "asd_analysis", "savepoint", "clusters.RData"))






# find NN for each cell:
library( RcppAnnoy )
featureMatrix <- pca$x; k_nn <- 50
annoy <- new( AnnoyEuclidean, ncol(featureMatrix) )
for( i in 1:nrow(featureMatrix) ) 
  annoy$addItem( i-1, featureMatrix[i,] )
annoy$build( 50 ) # builds a forest  of n_trees trees. More trees gives higher precision when querying.
nn_cells <- t( sapply( 1:annoy$getNItems(), function(i) annoy$getNNsByItem( i-1, k_nn) + 1 ) )
nndists_cells <- sapply( 1:ncol(nn_cells), function(j) sqrt( rowSums( ( featureMatrix - featureMatrix[ nn_cells[,j], ] )^2 ) ) )
rm(featureMatrix, annoy)

# cluster on nearest neighbor graph (Louvain):
adj <- Matrix(0, nrow = nrow(pca$x), ncol = nrow(pca$x)) # has to be sparse, otherwise takes 80 GB of RAM
for(i in 1:ncol(nn_cells))
  adj[ cbind(1:nrow(pca$x), nn_cells[, i]) ] <- 1
for(i in 1:ncol(nn_cells))
  adj[ cbind(nn_cells[, i], 1:nrow(pca$x)) ] <- 1
cl_louvain <- cluster_louvain(  graph_from_adjacency_matrix(adj, mode = "undirected") )

# merge clusters that are separated by patient heterogeneity:
tmp_clusters <- cl_louvain$membership
tmp_clusters <- case_when(tmp_clusters %in% c(5, 6, 8, 1, 10, 20, 2) ~ 5, TRUE ~ tmp_clusters) # excit. Ns
tmp_clusters <- case_when(tmp_clusters %in% c(11, 15, 19) ~ 11, TRUE ~ tmp_clusters) # astrocytes
tmp_clusters <- case_when(tmp_clusters %in% c(3, 9, 18) ~ 3, TRUE ~ tmp_clusters) # OPCs





# Louvain clusters 
p_louv <- ggplot()+ coord_fixed() +
  geom_point(data = data.frame(umap_euc, cl=factor(tmp_clusters)),
             aes(X1, X2, col = cl), size = .1) +
  geom_label(data = group_by(data.frame(umap_euc, cl=factor(tmp_clusters)), cl) %>%summarise(X1=mean(X1), X2=mean(X2)), 
             aes(X1, X2, label = cl))
p_louv

# clusters from paper
p_paper <- ggplot()+ coord_fixed()+
  geom_point(data =data.frame(cell = colnames(counts), umap_euc) %>%
               left_join(select(cellinfo, cell, cluster), by="cell"),
             aes(X1, X2, col = cluster), size = .1) +
  geom_label(data = data.frame(cell = colnames(counts), umap_euc) %>%
               left_join(select(cellinfo, cell, cluster), by = "cell") %>% group_by(cluster) %>%
               summarise(X1=mean(X1), X2=mean(X2)),
             aes(X1, X2, label = cluster))
p_paper

# 
# save(list = c("cl_louvain", "tmp_clusters", "nn_cells", "nn_inothercluster"),
#      file = file.path("~", "asd_analysis", "savepoint", "clusters.RData"))






# Doublets and ambiguous cells ----------------------------------

# load (or re-execute everything in this section):
load(file.path("~", "asd_analysis", "savepoint", "doublets.RData"))







# number of NN from different cluster:
nn_inothercluster <- colSums(
  matrix(tmp_clusters[ t(nn_cells) ],
         ncol = nrow(nn_cells))   != 
  matrix(rep(tmp_clusters, each = ncol(nn_cells)),
         ncol = nrow(nn_cells)) )
 

# in silico doublets: randomly draw cells from different clusters and pool their UMIs to form a "synthetic" doublet:

cellsA <- sample(1:ncol(counts), 50000)
cellsB <- rep(NA, 50000)
smpA <- cellinfo$sample[cellsA]
clA  <- tmp_clusters[cellsA]
tmp <- data.frame(smpA, clA) %>% group_by(smpA, clA) %>% tally
for(i in 1:nrow(tmp)) {
  is_smp <- cellinfo$sample[cellsA] == tmp$smpA[i]  
  is_cl  <- tmp_clusters[cellsA] == tmp$clA[i]
  # sample amongst cells from same sample and different cluster:
  cellsB[ is_smp & is_cl ] <- base::sample(
    x = which(cellinfo$sample == tmp$smpA[i]   & !tmp_clusters == tmp$clA[i]),
    size = tmp$n[i],
    replace = T) # in case one cluster is larger than all others combined
}

doublet_raw <- Ccounts[, cellsA] + Ccounts[, cellsB]
doublet_pcs <- predict(pca,
                       newdata = sqrt( (t(doublet_raw) / colSums(doublet_raw))[, is_informative] ))


# Alternative 1 (clearer):
a <- FNN::get.knn(rbind(pca$x, doublet_pcs), k = 50)
nn_doublets <- a$nn.index
nndists_doublets <- a$nn.dist

# Alternative 2 (faster):
library( RcppAnnoy )
featureMatrix <- rbind(pca$x, doublet_pcs); k_nn <- 50
annoy <- new( AnnoyEuclidean, ncol(featureMatrix) )
for( i in 1:nrow(featureMatrix) ) 
  annoy$addItem( i-1, featureMatrix[i,] )
annoy$build( 50 ) # builds a forest  of n_trees trees. More trees gives higher precision when querying.
nn_doublets <- t( sapply( 1:annoy$getNItems(), function(i) annoy$getNNsByItem( i-1, k_nn) + 1 ) )
nndists_doublets <- sapply( 1:ncol(nn_doublets), function(j) sqrt( rowSums( ( featureMatrix - featureMatrix[ nn_doublets[,j], ] )^2 ) ) )
rm(featureMatrix, annoy)


# percentage of synthetic doublets in neighborhood for each cell:
dblts_perc <- rowMeans( nn_doublets > ncol(counts) )[ 1:ncol(counts) ]



# Run UMAP with Annoy's output
ump2 <- uwot::umap( NULL, nn_method = list( idx=nn_doublets, dist=nndists_doublets), 
                    n_threads=40, spread = 15, verbose=TRUE )

is_synth <- 1:nrow(ump2) > nrow(pca$x)






# save(list = c("nn_doublets", "nndists_doublets", "cellsA", "cellsB",
#                 "dblts_perc", "is_synth", "ump2"),
#      file = file.path("~", "asd_analysis", "savepoint", "doublets.RData"))








  



# Cluster contributions to doublets ---------------------------------------

gg <- ggplot_build(p_louv)
cl_cols <- unique(gg$data[[2]][c("label","colour")])

plot_grid(plotlist = c(list(p_louv),
lapply(c(21,1, 9, 6, 13), function(cl){
  ggplot()+coord_fixed()+
  geom_point(data=data.frame(ump2[!is_synth,]), aes(X1, X2), col="grey", size=.05)+
  geom_point(data = data.frame(ump2[is_synth,],
                               cl_contributed = tmp_clusters[cellsA] == cl | tmp_clusters[cellsB] == cl),
             aes(X1, X2, col = cl_contributed), size=.05) +
   scale_color_manual(values = c(`FALSE`="grey", `TRUE`=cl_cols[cl_cols$label==cl,"colour"])) +
    
    geom_label(data = data.frame(cell=colnames(counts), ump2[!is_synth,], cluster = tmp_clusters)%>%
                 filter(cluster == cl) %>% summarise(X1=mean(X1), X2=mean(X2)),
               aes(X1, X2, label = cl),
               fontface = "bold")+ theme(legend.position = "none") +
   scale_fill_manual(values = cl_cols[cl_cols$label==cl,"colour"]) +
  
   ggtitle(paste0("Doublets with contribution from cluster ", cl)) + theme(legend.position = "none")
})
 ))



# synthetic doublets have synthetic doublets in their neighborhood
a <- FNN::get.knn(rbind(pca$x, doublet_pcs), k = 50)
data.frame(perc_dbl = rowMeans( a$nn.index > nrow(pca$x) ), is_synth = 1:nrow(a$nn.index) > nrow(pca$x)) %>%
  ggplot() + geom_histogram(aes(perc_dbl, fill = is_synth), alpha=.3) + coord_cartesian(ylim = c(0, 14000))+
  facet_wrap(~is_synth)


plot(rowMeans( a$nn.index > nrow(pca$x) ), pch=20, cex=.4); abline(h = nrow(pca$x))






# Markers and doublets ----------------------------------------------------


# "MOC","PTPRC" 
markers <- c(astro="AQP4",
             astro="GFAP",
             oligod="PLP1", 
             schwann="MPZ",
             Tcell="SKAP1",
             OPC="TNR",
             OPC="PDGFRA",
             Endothel="VWF",
             pericytes ="PDGFRB",
             microglia="CD68",
             microglia="CD4",
             neurons="SYT1",
            stroma= "LAMA2",
            motorNeuron ="THY1",
            motorNeuron = "NEFL")

for(i in 1:length(markers)){
  g <- markers[i]
png(file.path("~", "asd_analysis","marker_umap", paste0(names(markers)[i], "_",g,".png")),
    width = 800, height = 800)
p <- data.frame(umap_euc, cellinfo, Gene = Tcounts[, g], sfs) %>%
  ggplot(aes(X1, X2, col = Gene / sfs / mean(1/sfs))) + geom_point(size=.1)+coord_fixed()+
  col_pwr_trans(1/2, g) + ggtitle(names(markers)[i])
print(p)
dev.off()
}






# sleepwalk: markers -----------------------------------------------------------------
# does not work in RStudio server yet.
library(sleepwalk)

getHighGenes <- function(marked){
  
  # If no genes are marked, clear output, do nothing else
  if( length(marked) == 0 ) {
    return( "" )
  }
  
  #for each gene we calculate mean and standart deviation in marked cells and all other cells
  df <- data.frame(
    meanMarked   =  apply( data$expr[ marked, ], 2, mean ),
    sdMarked     =  apply( data$expr[ marked, ], 2, sd ),
    meanUnmarked =  apply( data$expr[ -marked, ], 2, mean ),
    sdUnmarked   =  apply( data$expr[ -marked, ], 2, sd )
  )
  #separation score
  df$sepScore <- ( df$meanMarked - df$meanUnmarked ) / pmax( df$sdMarked + df$sdUnmarked, 0.002 )
  
  # round to two decimal places
  df <- round(df * 100)/100
  
  #print top genes  
  print(head( df[ order( df$sepScore, decreasing=TRUE ), ], 15 ))
  
  #make a plot of gene expression
  topGene <- rownames(df)[which.max(df$sepScore)]
  pl <- ggplot() + geom_point(aes(x = data$um[, 1], y = data$um[, 2], colour = data$expr[, topGene])) +
    labs(colour = topGene) + scale_color_gradient(low = "Yellow", high = "Red")
  
  print(pl)
}

sleepwalk(umap_euc, pca$x[, 1:20], pointSize = 2.5, on_selection = getHighGenes)











# smooth LFCs -------------------------------------------------------------
# For each cell, I look into it's neighborhood and compare ASD and control for
# a selected gene. Doing this for all cells gives me something like a smoothed
# logFoldChange that might be interesting to work with (use to cluster genes, etc.).




true_ids <- 1:ncol(counts)

# for every cell, find nearest control cells:
controlNNs <- get.knnx(
  data = pca$x[cellinfo$diagnosis == "Control",],
  query =pca$x,
  k = 50)
# indices should reference all cells, not just Controls:
controlNNs <- matrix(true_ids[cellinfo$diagnosis == "Control"][controlNNs$nn.index], ncol = 50)
# add index of the cell itself as first column:
controlNNs <- cbind(true_ids, controlNNs)

# same for neighboring ASD cells:
asdNNs <- get.knnx(
  data = pca$x[cellinfo$diagnosis == "ASD",],
  query =pca$x,
  k = 50)
asdNNs <- matrix(true_ids[cellinfo$diagnosis == "ASD"][asdNNs$nn.index], ncol = 50)
asdNNs <- cbind(true_ids, asdNNs)

# plot a cell (green) with its control neighbors (blue) and it's ASD neighbors (red):
i <- sample(1:ncol(counts), 1)
ggplot() + coord_fixed() + 
  geom_point(data= data.frame(umap_euc), aes(X1, X2), col="grey", size=.1) +
  geom_point(data= data.frame(umap_euc[controlNNs[i,-1],]), aes(X1, X2), col="blue", size=.5) +
  geom_point(data= data.frame(umap_euc[    asdNNs[i,-1],]), aes(X1, X2), col="red", size=.5) +
  geom_point(data= data.frame(umap_euc[  true_ids[i],, drop=FALSE]), aes(X1, X2), col="green", size=1) +
  ggtitle(paste0(cellinfo$diagnosis[i], " cell (green dot)")) +
  theme(plot.title = element_text(colour = c("ASD"="red", "Control"="blue")[cellinfo$diagnosis[i]]))



i <- 28247

sfs[controlNNs[i,]]
Tcounts[controlNNs[i,], "ZNF770"]

rbind(
  # cell itself
  data.frame(diagnosis = cellinfo$diagnosis[i], sfs=sfs[i], umi = Tcounts[i, "ZNF770"]),
  # neighbors:
  data.frame(
    diagnosis = rep(c("Control", "ASD"), each=50),
    sfs = c(sfs[controlNNs[i,-1]], sfs[asdNNs[i,-1]]),
    umi = c(Tcounts[controlNNs[i,-1], "ZNF770"], Tcounts[asdNNs[i,-1], "ZNF770"]))
  ) %>%
  ggplot(aes(sfs, umi, col = diagnosis))+geom_point()
 

 






# DESeq -------------------------------------------------------------------
library(DESeq2)
library(BiocParallel)

# visualize dirty cells we clean away:
tmp <- data.frame(umap_euc,
                  diagnosis = cellinfo$diagnosis,
                  clean = dblts_perc < 3/50  & nn_inothercluster < 1,
                  Gene = Tcounts[, "TTF2"] / sfs/mean(1/sfs),
                  cl = factor(tmp_clusters))
ggplot() + coord_fixed()+
  geom_point(data=filter(tmp,  clean), aes(X1, X2, col = cl), size=.1) +
  geom_point(data=filter(tmp, !clean), aes(X1, X2), col = "black", size=.1) +
  geom_label(data=group_by(tmp, cl) %>% summarise(X1=mean(X1), X2=mean(X2)), aes(X1, X2, label=cl))

tmp <- as.matrix(table(sample=cellinfo$sample, clean = dblts_perc < 3/50  & nn_inothercluster < 1))
data.frame(sample = rownames(tmp), dirtyProportion = tmp[,1] / (tmp[,1] + tmp[,2])) %>% left_join(sampleTable, by="sample") %>% ggplot(aes(sample, dirtyProportion, col = diagnosis))+geom_point()



# compute for a single cluster
sel <-Tcounts[, "SYT1"] > 1 & Tcounts[, "CUX2"] > 0  & dblts_perc < 3/50  & nn_inothercluster < 1
pseudobulks <- as.matrix(t( fac2sparse(cellinfo$sample[sel]) %*% t(Ccounts[, sel]) ))
coldat <- filter(sampleTable, sample %in% colnames(pseudobulks)) %>% 
  mutate(individual = factor(individual),
         diagnosis = factor(diagnosis, levels = c("Control", "ASD")),
         region    = factor(region))
rownames(coldat) <- coldat$sample

dds <- DESeqDataSetFromMatrix( pseudobulks,
                               coldat[colnames(pseudobulks), ],
                               design = ~ sex + region + age + diagnosis )
# For cluster 5, I tested that we do not need interactions between sex, region and diagnosis. I used
# DESeq's LTR for this (see mail to Simon at mid-September 2019).
dds <- DESeq(dds, 
             parallel=TRUE, BPPARAM=MulticoreParam(20))
res_df <- results(dds, name = "diagnosis_ASD_vs_Control") %>% as.data.frame() %>% rownames_to_column("Gene")
table(res_df$padj < .1)  
res_df %>% arrange(padj) %>% head(n=20)



cors <- cor(
 sqrt( Tcounts[sel, "EPB41L5"] / sfs[sel] ),
 as.matrix(sqrt( Tcounts[sel, ] / sfs[sel] ))
)




# correlation-like stuff --------------------------------------------------
clean <- dblts_perc < 3/50  & nn_inothercluster < 1

s_50 <- rowSums( matrix(sfs[ nn_cells ], ncol = 50) )   
syt1_50 <- rowSums( matrix(Tcounts[, "SYT1"][ nn_cells ], ncol = 50) )   / 50
cux2_50 <- rowSums( matrix(Tcounts[, "CUX2"][ nn_cells ], ncol = 50) )   / 50
ttf2_50 <- rowSums( matrix(Tcounts[, "TTF2"][ nn_cells ], ncol = 50) )   / 50

data.frame(umap_euc, s_50,
           syt1_50,
           cux2_50,
           ttf2_50,
           clean,
           diagnosis = cellinfo$diagnosis,
           syt1_raw = Tcounts[, "SYT1"],
           ttf2_raw = Tcounts[, "TTF2"]
           ) %>%
  gather(Gene, knn, syt1_50, cux2_50, ttf2_raw) %>% 
  ggplot(aes(X1, X2, col= knn / s_50 / mean(1/s_50)))+coord_fixed()+
  geom_point(size=.1) + facet_wrap(~ diagnosis + Gene) +
  scale_color_gradientn(
    trans = power_trans(1/2),
    colours = rev(rje::cubeHelix(100))[5:100],
    na.value = adjustcolor("grey", alpha.f = .4),
    labels = semi_scientific_formatting
  )




# markers for cortical layers, from this paper (Fig. 6):
#   Large-Scale Cellular-Resolution Gene Profiling in Human Neocortex Reveals
#   Species-Specific Molecular Signatures
#   Zeng, Shen, ..., Kleinman, Jones
#   Cell 2012
l1 <- c("NDNF", # aka "C4orf31"
        "CHRNA7_ENSG00000175344", "CHRNA7_ENSG00000274542", # CHRNA7 exists twice
        "CNR1","CXCL14","RELN","INPP4B")
l23<-c("LAMP5", # aka: "C20orf103",
  "GSG1L","IGSF11","KCNIP2","PVRL3","RASGRF2","SYT17","WFS1","C1QL2",
  "CARTPT","CALB1","CUX2","ATP2B4","CBLN2","CCK","FXYD6","PENK","CACNA1E","KCNH4",
  "SCN3B","COL24A1","CRYM","TPBG","BEND5","COL6A1","PRSS12","SCN4B","SYT2","LGALS1",
  "MFGE8","SV2C","SNCG")

l4 <- c("RORB","CACNG5","CHRNA3","GRIK4","KCNIP1","PDYN")

l5 <- c(
  "TRIB2","CPNE7","ETV1","FAM3C","TOX","VAT1L","KIAA1456","HTR2C"
)

l56 <- c("PCDH20_ENSG00000280165", "PCDH20_ENSG00000197991", # PCDH20 exists twice in our data
  "B3GALT2","KCNK2","PCP4","PDE1A","RPRM","RXFP1","GABRA5","KCNA1"
)

l6 <- c(
  "CDH24","CYR61","FOXP2","NTNG2","SYT10","SYT6","TH","TLE4","TMEM163","AKR1C2","AKR1C3","ANXA1","NPY2R","OPRK1","PCDH17","SEMA3C","SYNPR"
)

l6b_wm <- c( # so called "subplate neurons". Interstitial neurons are from white matter (wm)
  "ADRA2A","CTGF","NR4A2"
)


the_markers <- c(l1, l23, l4, l5, l56, l6, l6b_wm)

# from Schirmer Science paper:
interneuron_markers <- c("GAD1", "GAD2")

gene_means <- rowMeans( norm_counts[, tmp_clusters == 5] )
gene_vars <- rowVars_spm( norm_counts[, tmp_clusters == 5] )
frequent <- rowSums(norm_counts[, tmp_clusters == 5] != 0)  >  .01 * sum(tmp_clusters == 5)
tmp <- data.frame(gene = names(gene_means), gene_means, gene_vars,
                  is_marker = names(gene_means) %in% the_markers,
                  frequent, stringsAsFactors = F)  %>%
 mutate(use_marker = is_marker & frequent & gene_vars/gene_means > 1.5 * mean(1/sfs[tmp_clusters==5]))
ggplot() +
 geom_point(data = filter(tmp, !is_marker), aes(gene_means, gene_vars/gene_means), size=.1, col="grey")+
 geom_point(data = filter(tmp, is_marker & !use_marker), aes(gene_means, gene_vars/gene_means), size=.5, col="black")+
 geom_point(data = filter(tmp, use_marker), aes(gene_means, gene_vars/gene_means), size=.5, col="red")+
  scale_x_log10()+scale_y_log10() +
  geom_hline(yintercept = mean(1/sfs[tmp_clusters==5]))+
  geom_hline(yintercept = 1.5 * mean(1/sfs[tmp_clusters==5]), linetype="dashed", col="red")

# Correlation:
markers_use <- the_markers[the_markers %in% (filter(tmp, use_marker) %>% pull(gene))]
cors <- cor(as.matrix(t(sqrt(norm_counts[markers_use,
              tmp_clusters == 5 & clean] ))))
diag(cors) <- NA
library(RColorBrewer)
pheatmap::pheatmap(cors,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                             "RdBu")))(100),
                   cluster_rows = F, cluster_cols = F)


# where are the Layer4 Neurons?
data.frame(umap_euc,
           clean,
           sfs,
           diagnosis = cellinfo$diagnosis,
           louv = tmp_clusters,
           as.matrix(Tcounts[, l4[l4 %in% markers_use]])
           ) %>% 
  gather(Gene, UMI, -X1, -X2, -clean, -sfs, -diagnosis, - louv) %>%
  ggplot(aes(X1, X2, col = UMI/sfs/mean(1/sfs)))+geom_point(size=.1) + coord_fixed()+
  col_pwr_trans(1/10)+
  facet_wrap(~Gene)

# they use RORB as marker for L4, but it's correlated highly with these l5/6 markers:
l56_and_rorb <- c("RORB","TOX","KIAA1456","PDE1A","RXFP1","FOXP2")
  
data.frame(umap_euc,
           sfs,
           as.matrix(Tcounts[, l56_and_rorb])) %>%
  gather(Gene, UMI, -X1, -X2, -sfs) %>%
  ggplot(aes(X1, X2, col = UMI/sfs/mean(1/sfs)))+geom_point(size=.1) + coord_fixed()+
  col_pwr_trans(1/10)+
  facet_wrap(~Gene)






# more DESeq stuff --------------------------------------------------------


# compare ASD vs control for several clusters:
res_clean <- lapply(c(5, 16, 21, 4, 11, 13, 3, 14, 7, 23), function(cl){
  print(cl)
  sel <- tmp_clusters== cl      & dblts_perc < 3/50  & nn_inothercluster < 1
  
  pseudobulks <- as.matrix(t( fac2sparse(cellinfo$sample[sel]) %*% t(Ccounts[, sel]) ))
  coldat <- filter(sampleTable, sample %in% colnames(pseudobulks)) %>% 
    mutate(individual = factor(individual),
           sex = factor(sex),
           diagnosis = factor(diagnosis, levels = c("Control", "ASD")),
           region    = factor(region))
  rownames(coldat) <- coldat$sample
  dds <- DESeqDataSetFromMatrix( pseudobulks, coldat[colnames(pseudobulks), ],
                                 design = ~ sex + age + region + diagnosis )
  dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(20))
  res_df <- results(dds, name = "diagnosis_ASD_vs_Control") %>% as.data.frame() %>% rownames_to_column("Gene") 
  list(cluster = cl, ncells = sum(sel), res = res_df ) 
})

names(res_clean) <- unlist(lapply(res_clean, function(x) x$cluster))


plot_grid(plotlist = lapply(names(res_clean), function(cl){
 data.frame(padj_clean=res_clean[[cl]]$res$padj,
           padj_dirty=res_dirty[[cl]]$res$padj) %>% ggplot(aes(-log10(padj_dirty), -log10(padj_clean)))+
  geom_point(size=.1)+coord_fixed() + geom_abline() + geom_vline(xintercept = 1, lty=2, col="red")+
  geom_hline(yintercept = 1, lty=2, col="red")+ ggtitle(cl)
 
}) )


# Compare to sfari database
sfari <- read_csv(file.path("~", "asd_analysis",
                            "SFARI-Gene_genes_08-29-2019release_09-24-2019export.csv")) %>%
  rename_all(make.names)
in_database <- sfari$gene.symbol[ sfari$gene.symbol %in% gene_info$V2 ]


tmp <- lapply(names(res_clean), function(cl){
  print(cl)
  degs <-    res_dirty[[cl]]$res %>% filter(padj < .1) %>% pull(Gene)
  in_test <- res_dirty[[cl]]$res %>% filter(!is.na(padj)) %>% pull(Gene)
  if(length(degs)==0){p_dirty <- NA}else{
  p_dirty <- fisher.test(matrix(table( in_test %in% degs, in_test %in% in_database),
                     ncol=2,
                     dimnames = list(is_deg = c("no","yes"), in_database = c("no","yes"))))$p.value
  }
  degs <-    res_clean[[cl]]$res %>% filter(padj < .1) %>% pull(Gene)
  in_test <- res_clean[[cl]]$res %>% filter(!is.na(padj)) %>% pull(Gene)
  if(length(degs)==0){p_clean <- NA}else{
  p_clean <- fisher.test(matrix(table( in_test %in% degs, in_test %in% in_database),
                     ncol=2,
                     dimnames = list(is_deg = c("no","yes"), in_database = c("no","yes"))))$p.value
  }
  return(data.frame(cluster=cl, p_dirty = p_dirty, p_clean = p_clean, stringsAsFactors = F))
}) %>% bind_rows

tmp %>% ggplot(aes(-log10(p_dirty), -log10(p_clean)))+geom_point()+geom_abline() + ggtitle("DEG enrichment")

data.frame(
  cluster = names(res_clean),
  ncells_clean = lapply(res_clean, function(x) x$ncells) %>% unlist,
  ncells_dirty = lapply(res_dirty, function(x) x$ncells) %>% unlist
) %>% left_join(tmp) %>% head









res_df %>% filter(padj < .1, baseMean > 50) %>% arrange(desc(abs(log2FoldChange))) %>% head(n=20)

Gene   baseMean log2FoldChange     lfcSE      stat       pvalue        padj
1  MTND2P28   60.25608      2.7302927 0.6788128  4.022159 5.766718e-05 0.012415743
2     HSPB1   53.25841     -1.8056395 0.5872707 -3.074629 2.107647e-03 0.060330707
3    MT-ND3  870.50560      1.4187626 0.3279703  4.325887 1.519194e-05 0.006305464
4   MT-ND4L  167.39401      1.3428971 0.3245225  4.138071 3.502378e-05 0.009425774
5    MT-ND4 1653.17914      1.3137104 0.2869792  4.577720 4.700708e-06 0.003935798
6  MTATP6P1  361.22570      1.0486883 0.3376722  3.105640 1.898675e-03 0.057480286
7    MT-CO3 2461.32824      1.0356872 0.2742279  3.776739 1.588948e-04 0.020123562
8    MT-CO2 1831.71986      0.9955536 0.2746291  3.625084 2.888675e-04 0.026873597
9    MT-ND1  658.56180      0.9821369 0.2345878  4.186650 2.831022e-05 0.008585289
10  MT-ATP6  874.24854      0.9531860 0.2803692  3.399752 6.744688e-04 0.038315943
11   MT-ND2  800.02263      0.9369508 0.2443916  3.833810 1.261736e-04 0.017607065
12     TTF2   50.01172     -0.8898189 0.1978098 -4.498355 6.848117e-06 0.004549986
13   MT-CYB  718.85866      0.8116936 0.2152989  3.770077 1.631970e-04 0.020272786
14    ZMYM3   57.76922      0.8077571 0.2202002  3.668286 2.441821e-04 0.024838436
15      KIT  169.95008      0.7833860 0.1813739  4.319178 1.566115e-05 0.006305464
16    VMA21   68.14624      0.7507723 0.2401887  3.125760 1.773462e-03 0.055716836
17    CLIP3  143.12863      0.7417990 0.1567059  4.733701 2.204624e-06 0.003524178
18     PLK2   55.06825      0.7178509 0.1663288  4.315856 1.589859e-05 0.006305464
19    CXXC4   52.21341      0.7035908 0.1437565  4.894323 9.864475e-07 0.002973350
20    KLHL9   61.72171      0.6984751 0.2225714  3.138207 1.699847e-03 0.054623451






# Plot individual genes
g <- "MT-ND3"
plotCounts(dds, g, intgroup = c("sex", "region", "diagnosis"))
data.frame(umap_euc, cellinfo, Gene = Tcounts[, g], sfs, sel) %>%
   # filter(sel) %>%
  ggplot(aes(X1, X2, col = Gene / sfs / mean(1/sfs))) + geom_point(size=.1)+coord_fixed()+
  col_pwr_trans(1/2, g) + facet_wrap(~ region + diagnosis)









# Investigate gene-gene correlations (maybe useful to see if subpopulations of cells exist or not):
deg_cors <- cor(  as.matrix( t(sqrt(norm_counts[degs, sel])) ) )
hist(deg_cors, 100)
# adjacency matrix: 
deg_adj <- 0 + (deg_cors > .2)

neighborless <- rowSums( deg_cors > .2 ) <= 1
deg_cl <- cluster_louvain( graph_from_adjacency_matrix(deg_adj[!neighborless, !neighborless], mode = "undirected") )

deg_umap <- uwot::umap( 1-deg_cors[!neighborless, !neighborless], spread = 10, n_threads = 10)
data.frame(deg_umap, cl = factor(deg_cl$membership)) %>% ggplot(aes(X1, X2, col=cl))+geom_point()


groups(deg_cl)$`1`  # investigate further?








# NRGN neurons ------------------------------------------------------------



sel <- cellinfo$cluster == "Neu-NRGN-II"  #grepl("NRGN", cellinfo$cluster)
pseudobulks <- as.matrix(t( fac2sparse(cellinfo$sample[sel]) %*% t(Ccounts[, sel]) ))
coldat <- filter(sampleTable, sample %in% colnames(pseudobulks)) %>% 
  mutate(individual = factor(individual),
         diagnosis = factor(diagnosis, levels = c("Control", "ASD")),
         sex = factor(sex),
         region    = factor(region))
rownames(coldat) <- coldat$sample

dds <- DESeqDataSetFromMatrix( pseudobulks,
                               coldat[colnames(pseudobulks), ],
                               design = ~ sex + region + age+ diagnosis )
dds <- DESeq(dds, 
             parallel=TRUE, BPPARAM=MulticoreParam(20))
res_NRGN2 <- results(dds, name = "diagnosis_ASD_vs_Control") %>% as.data.frame() %>% rownames_to_column("Gene")



data.frame(
  Gene =     rownames(res_NRGN1),
  # p_pooled = res_pooledNRGN$padj,
  p_nrgn1  = res_NRGN1$padj,
  p_nrgn2  = res_NRGN2$padj
) %>%
  ggplot(aes(-log10(p_nrgn1), -log10(p_nrgn2)))+geom_point()+coord_fixed() + geom_abline()


# ASD/cntrl and clusterI/II are perfectly mixed, II is just 2.5x larger than I:
filter(cellinfo, grepl("NRGN", cluster)) %>% select(cluster, diagnosis) %>% table()








# DE playground -----------------------------------------------------------

# DE testing has to follow simple assumptions to be feasible, as there are
# infinitely many possible distributions the counts could have in controls and treatments,
# respectively.
# For example, perhaps in my 20 control patients a gene is distributed according to
# a NB in some and according to a composite of two NBs in other samples, and another
# gene according to a composite of three NBs with different means.
# It is infeasible to model this with NB/composite-NB fits.
# Instead, let's think about genes we would be interested in as DEGs in scRNAseq:
# such a gene would in the treated samples follow distributions (NB / composite-NBs)
# whose parameters are sampled from a different entity - e.g. they could be 
# composite-NBs composed of three instead of two NBs, or they could be
# composite-NBs around the same two means, but more cells with the higher mean, etc..
# Again, it is completely infeasible to model this or do fits with EM algorithm or whatever,
# instead we have to come up with simple assumptions that describe this well.

# 


g <- "RGS4"

data.frame(
  sample = aggregate(x = sfs[sel], by = list(sample=cellinfo$sample[sel]), FUN = median)[, 1],
  sf   = aggregate(x = sfs[sel], by = list(sample=cellinfo$sample[sel]), FUN = median)$x,
  mean = aggregate(x = norm_counts[g, sel], by = list(sample=cellinfo$sample[sel]), FUN = mean)$x,
  sdev   = aggregate(x = norm_counts[g, sel], by = list(sample=cellinfo$sample[sel]), FUN = sd)$x,
  stringsAsFactors = F
) %>% left_join( select(sampleTable, sample, diagnosis), by = "sample" ) %>% 
  ggplot() + geom_point(aes(mean, sdev, col = diagnosis))  + scale_x_log10()  + scale_y_log10()






