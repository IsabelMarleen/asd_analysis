# This script crashes with 16 GB or less of RAM.

library( tidyverse )
library( Matrix )
library( irlba )
library( uwot )
library( cowplot )
scr_dir <- "/home/frauhammer/sc_methods_dev/src/"
source(file.path(scr_dir, "functions_universal.R"))


path <- "/home/frauhammer/sds_copy/ASD/"
cellinfo <- read.delim( file.path( path, "meta.txt" ), stringsAsFactors=FALSE )
counts <- readMM( file.path( path, "matrix.mtx" ) )
rownames(counts) <- read.delim( file.path( path, "genes.tsv" ), header=FALSE )$V2
colnames(counts) <- readLines( file.path( path, "barcodes.tsv" ) )


sampleTable <-
  cellinfo %>% select( sample : RNA.Integrity.Number ) %>% unique
sampleTable

# extracting gene expression is much faster in column-sparse format:
Tcounts <- as(t(counts), "dgCMatrix")   #  fast:   Tcounts[, "SYN1"]


# find informative genes
sfs <- colSums(counts)
norm_counts <- t(t(counts) / colSums(counts))
poisson_vmr <- mean(1/sfs)
gene_means <- rowMeans( norm_counts )
gene_vars <- rowVars_spm( norm_counts )
cells_expressing <- rowSums( counts != 0 )
is_informative <- gene_vars/gene_means > 1.5 * poisson_vmr  &  cells_expressing > 100

plot(gene_means, gene_vars/gene_means, pch=".", log = "xy")
points(gene_means[is_informative], (gene_vars/gene_means)[is_informative], pch=".", col = "red" )



# PCA and UMAP
pca <- irlba::prcomp_irlba( x = sqrt(t(norm_counts[is_informative,])),
                            n = 50,
                            scale. = TRUE)
umap_euc <- uwot::umap( pca$x, spread = 2)

umap_cos <- uwot::umap( pca$x, metric = "cosine", spread = 2)



# Clusters from paper -----------------------------------------------------

df <- data.frame(cell = colnames(norm_counts), umap_euc, stringsAsFactors = F) %>%
  left_join(cellinfo, by = "cell") %>%
  mutate(individual = factor(individual))
# some clusters from the paper are mixing in UMAP, i.e. are not robust to preprocessing, so we merge them:
df <- df %>% mutate(cluster = case_when(
  grepl("^Neu-NRGN-", cluster) ~ "Neu_NRGN",
  grepl("^AST-", cluster) ~ "AST",
  TRUE ~ cluster
))

labels_df <- df %>% group_by(cluster) %>% summarise(X1 = mean(X1), X2=mean(X2))

p_cl <-   ggplot() + 
  geom_point(data=df,
             aes(X1, X2, col = cluster),
             size = .1) + coord_fixed()+
  geom_label(data = labels_df,
             aes(X1, X2, col = cluster, label = cluster),
             fontface = "bold")+ theme(legend.position = "none")
p_cl





# Clusters and doublets ---------------------------------------------------


g <- "SYT1"    # Neuron marker
g <- "SLC1A2"  # Astrocyte marker

ggplot()+
  geom_point(data = data.frame(umap_euc, Gene = Tcounts[, "SLC1A2"], sfs = sfs),
             aes(X1, X2, col = Gene / sfs / mean(1/sfs)), size=.1) +
  coord_fixed() + col_pwr_trans(1/2) +
  geom_label(data = labels_df,
             aes(X1, X2, label = cluster),
             fontface = "bold")



