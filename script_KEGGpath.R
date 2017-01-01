setwd('/home/Yulong/RESEARCH/BaiLinhan_array/Process/')

library('KEGGAPI')
library('clusterProfiler')


load('mergeGeneAnno.RData')
DEGMat <- read.csv('DEGMat.csv', stringsAsFactor = FALSE, row.names = 1)
mergeGeneAnno[, 3] <- sapply(strsplit(as.character(mergeGeneAnno[, 3]), split = ':', fixed = TRUE), '[[', 2)

## DEG genes
## threshold is P.Val < 0.05
thres <- 0.05
DEGProbes <- rownames(DEGMat)[DEGMat[, 'P.Value'] < 0.05]
DEGIDs <- mergeGeneAnno[mergeGeneAnno[, 'ProbeName'] %in% DEGProbes, 'KEGGID']

kk <- enrichKEGG(gene = DEGIDs, organism = 'mgr', pvalueCutoff = 1, qvalueCutoff = 1, universe = mergeGeneAnno[, 3], minGSSize = 5)
showNum <- 10
pdf('KEGGbar.pdf')
barplot(kk, showCategory = showNum)
dev.off()
pdf('KEGGdot.pdf')
dotplot(kk, showCategory = showNum)
dev.off()

## convert to probe IDs
kk <- as.data.frame(kk)
kk[, 8] <- sapply(kk[, 8], function(x) {
  splitVec <- unlist(strsplit(x, split = '/', fixed = TRUE))
  convVec <- mergeGeneAnno[mergeGeneAnno[, 'KEGGID'] %in% splitVec, 'ProbeName']
  convVec <- paste(convVec, collapse = '/')
  return(convVec)
})

write.csv(kk, 'KEGGenrich.csv')

