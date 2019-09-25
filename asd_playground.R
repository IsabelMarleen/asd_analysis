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
p_louv <- ggplot()+
  geom_point(data = data.frame(umap_euc, cl=factor(tmp_clusters)),
             aes(X1, X2, col = cl), size = .1) +
  geom_label(data = group_by(data.frame(umap_euc, cl=factor(tmp_clusters)), cl) %>%summarise(X1=mean(X1), X2=mean(X2)), 
             aes(X1, X2, label = cl))
p_louv

# clusters from paper
ggplot()+
  geom_point(data =data.frame(cell = colnames(counts), umap_euc) %>%
               left_join(select(cellinfo, cell, cluster), by="cell"),
             aes(X1, X2, col = cluster), size = .1) +
  geom_label(data = data.frame(cell = colnames(counts), umap_euc) %>%
               left_join(select(cellinfo, cell, cluster), by = "cell") %>% group_by(cluster) %>%
               summarise(X1=mean(X1), X2=mean(X2)),
             aes(X1, X2, label = cluster))


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
dblts_perc <- rowMeans( nn_doublets > nrow(pca$x) )[ 1:nrow(pca$x) ]



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








# DESeq -------------------------------------------------------------------
library(DESeq2)
library(BiocParallel)

# visualize dirty cells we clean away:
tmp <- data.frame(umap_euc, clean = dblts_perc == 0  & nn_inothercluster <= 1)
ggplot() + geom_point(data=filter(tmp, clean), aes(X1, X2), col = "grey", size=.1) + geom_point(data=filter(tmp, !clean), aes(X1, X2), col = "black", size=.1)



# compare ASD vs control for one cluster:
sel <- tmp_clusters==5  & dblts_perc == 0  & nn_inothercluster <= 1
pseudobulks <- as.matrix(t( fac2sparse(cellinfo$sample[sel]) %*% t(Ccounts[, sel]) ))
coldat <- filter(sampleTable, sample %in% colnames(pseudobulks)) %>% 
  mutate(individual = factor(individual),
         age = factor(age),
         diagnosis = factor(diagnosis, levels = c("Control", "ASD")),
         region    = factor(region))
rownames(coldat) <- coldat$sample

dds <- DESeqDataSetFromMatrix( pseudobulks,
                               coldat[colnames(pseudobulks), ],
                               design = ~ sex + age + diagnosis )
# For cluster 5, I tested that we do not need interactions between sex, region and diagnosis. I used
# DESeq's LTR for this (see mail to Simon at mid-September 2019).
dds <- DESeq(dds, 
             parallel=TRUE, BPPARAM=MulticoreParam(20))




# Plot individual genes
g <- "RGS4"
plotCounts(dds, g, intgroup = c("sex", "region", "diagnosis"))
data.frame(umap_euc, cellinfo, Gene = Tcounts[, g], sfs, sel) %>%
   filter(sel) %>%
  ggplot(aes(X1, X2, col = Gene / sfs / mean(1/sfs))) + geom_point(size=.1)+coord_fixed()+
  col_pwr_trans(1/2, g) + facet_wrap(~ region + diagnosis)





# Compare to sfari database

sfari <- read_csv(file.path("~", "asd_analysis",
                            "SFARI-Gene_genes_08-29-2019release_09-24-2019export.csv")) %>%
  rename_all(make.names)
in_database <- sfari$gene.symbol[ sfari$gene.symbol %in% gene_info$V2 ]

dds <- dds_32k
degs <- results(dds, name = "diagnosis_ASD_vs_Control") %>% as.data.frame() %>% rownames_to_column("Gene") %>%
  filter(padj < .1) %>% arrange(desc(log2FoldChange)) %>% pull(Gene)
in_test     <- results(dds, name = "diagnosis_ASD_vs_Control") %>% as.data.frame() %>%
                   rownames_to_column("Gene") %>% filter(!is.na(padj)) %>% pull(Gene) 

for_fisher <- matrix(table( in_test %in% degs, in_test %in% in_database), ncol=2,
       dimnames = list(is_deg = c("no","yes"), in_database = c("no","yes")))

fisher.test(for_fisher)




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


