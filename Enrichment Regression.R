er <- log2(
  table( cellTable$newclusters, cellTable$sample ) / nrow(cellTable) /
    outer( table( cellTable$newclusters ) / nrow(cellTable),
           table( cellTable$sample ) / nrow(cellTable) ) )
er[ er < -5] <- -5
er[ er > 5] <- 5
pheatmap::pheatmap( er[ , order(sampleTable$diagnosis) ], cluster_cols = FALSE )

cellTable %>%
  group_by( sample ) %>%
  summarise( ncells = n(), nAST = sum( newclusters == 2 ) ) %>%
  left_join( sampleTable ) %>%
ggplot +
  geom_point(aes( x=ncells, y= nAST/ncells, col=diagnosis))
