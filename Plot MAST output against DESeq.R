#Plotting q values from MAST output against padj from DESeq output
res[ order(res$padj), ]
results_lucas <- read.delim(file="~/Desktop/ASD/S4.csv", sep=";", dec=",")
results_lucas <- results_lucas[results_lucas$X...Cell.type == "L2/3", ]
res_L23 <- res[rownames(res) %in% results_lucas$Gene.name, ]

results1 <- left_join(as.tibble(results_lucas), 
          as.tibble(rownames_to_column(as.data.frame(res_L23), var="Gene.name")), by="Gene.name")
plot(log(results1$q.value), log(results1$padj), xlab="q.value from Schirmer", ylab="padj from DESeq")
abline(0,1)

res_HS6ST3 <- res2[rownames(res2) %in% results_lucas$Gene.name, ]

results <- left_join(as_tibble(results_lucas), 
                     as_tibble(res_HS6ST3, rownames="Gene.name"))
plot(log(results$q.value), log(results$padj2), xlab="q.value from Schirmer", ylab="padj from DESeq clustered 
     by HST...")
abline(0,1)


tibble(counts["HS6ST3",])
as_tibble(cellinfo)%>%
  add_column(count=counts["HS6ST3",])%>%
  group_by(sample) %>%
  summarise(auc=pROC::auc(pROC::roc(I(cluster=="L2/3"), count/UMIs)))

res[ order(res$padj), ]
res2[order(res$padj), ]
colnames(res2) <- c("baseMean2", "log2FoldChange2", "lfcSE2", "stat2", "pvalue2", "padj2")
results_joined <- left_join(as_tibble(res2, rownames="Gene.name"), 
                            as_tibble(res, rownames="Gene.name"))
plot(log10(results_joined$padj), log10(results_joined$padj2), xlab="Sorting by Paper", 
     ylab="Sorting by HS6ST3")
abline(0,1)

all(rownames(res[order(res$padj), ]) %in% rownames(res2[order(res$padj), ]))
