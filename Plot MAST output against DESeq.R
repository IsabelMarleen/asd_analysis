#Plotting q values from MAST output against padj from DESeq output

res[ order(res$padj), ]
results_lucas <- read.delim(file="~/Desktop/ASD/S4.csv", sep=";", dec=",")
results_lucas <- results_lucas[results_lucas$X...Cell.type == "L2/3", ]
res_L23 <- res[rownames(res) %in% results_lucas$Gene.name, ]
results_joined <- left_join(as.tibble(results_lucas), 
          as.tibble(rownames_to_column(as.data.frame(res_L23), var="Gene.name")), by="Gene.name")
plot(log(results_joined$q.value), log(results_joined$padj), xlab="q.value from Schirmer", ylab="padj from DESeq")



tibble(counts["HS6ST3",])
as_tibble(cellinfo)%>%
  add_column(count=counts["HS6ST3",])%>%
  group_by(sample) %>%
  summarise(auc=pROC::auc(pROC::roc(I(cluster=="L2/3"), count/UMIs)))

res[ order(res$padj), ]
res2[order(res$padj), ]
colnames(res2) <- c("baseMean2", "log2FoldChange2", "lfcSE2", "stat2", "pvalue2", "padj2")
results_joined <- left_join(as_tibble(rownames_to_column(as.data.frame(res2))), 
                            as_tibble(rownames_to_column(as.data.frame(res))))
plot(log10(results_joined$padj), log10(results_joined$padj2), xlab="Sorting by Paper", ylab="Sorting by HS6ST3")
abline(0,1)

all(rownames(res[order(res$padj), ]) %in% rownames(res2[order(res$padj), ]))
