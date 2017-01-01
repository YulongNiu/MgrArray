########################convert GO##################################
setwd('/home/Yulong/RESEARCH/BaiLinhan_array/Process/')

library('clusterProfiler')
library('stringr')
library('GO.db')

load('mergeProAnno.RData')
mergeProAnno[, 1] <- sapply(strsplit(as.character(mergeProAnno[, 1]), split = ':', fixed = TRUE), '[[', 2)
DEGMat <- read.csv('DEGMat.csv', stringsAsFactor = FALSE, row.names = 1)

rawGO <- read.csv('gene_association.PAMGO_Mgrisea', sep = '\t', skip = 24, stringsAsFactor = FALSE, header = FALSE)
rawGO <- rawGO[, c(11, 5, 9)]

## 1. remove NA character string
proGO <- rawGO[rawGO[, 1] != '', ]
## 2. keep MGG_**
mggNames <- sapply(strsplit(proGO[, 1], split = '|', fixed = TRUE), function(x) {
  x <- unlist(str_extract_all(x, 'MGG_.*'))
  return(x)
})
proGO <- cbind(unlist(mggNames),
               rep(proGO[, 2], sapply(mggNames, length)),
               rep(proGO[, 3], sapply(mggNames, length)))
## 3. add term annotation
termAnno <- sapply(proGO[, 2], function(x) {
  return(Term(GOTERM[[x]]))
})
proGO <- cbind(proGO, termAnno)
rownames(proGO) <- NULL


## DEG genes
## threshold is P.Val < 0.05
thres <- 0.05
DEGProbes <- rownames(DEGMat)[DEGMat[, 'P.Value'] < 0.05]
DEGIDs <- mergeProAnno[mergeProAnno[, 'ProbeName'] %in% DEGProbes, 'KEGGID']

##~~~~~~~~~~~~~~~~~~~~~~~~~BP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
proGOBP <- proGO[proGO[, 3] == 'P', ]
bp <- enricher(DEGIDs, pvalueCutoff = 1, universe = unique(proGO[, 1]), minGSSize = 5, qvalueCutoff = 1, TERM2GENE = data.frame(proGOBP[, 2:1]), TERM2NAME = data.frame(proGOBP[, c(2, 4)]))

showNum <- 10
pdf('BPbar.pdf', width = 10)
barplot(bp, showCategory = showNum)
dev.off()
pdf('BPdot.pdf', width = 10)
dotplot(bp, showCategory = showNum)
dev.off()

## convert to probe IDs
bp <- as.data.frame(bp)
bp[, 8] <- sapply(bp[, 8], function(x) {
  splitVec <- unlist(strsplit(x, split = '/', fixed = TRUE))
  convVec <- mergeProAnno[mergeProAnno[, 'KEGGID'] %in% splitVec, 'ProbeName']
  convVec <- paste(convVec, collapse = '/')
  return(convVec)
})

write.csv(bp, file = 'BPenrich.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~MF~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
proGOMF <- proGO[proGO[, 3] == 'F', ]
mf <- enricher(DEGIDs, pvalueCutoff = 1, universe = unique(proGO[, 1]), minGSSize = 5, qvalueCutoff = 1, TERM2GENE = data.frame(proGOMF[, 2:1]), TERM2NAME = data.frame(proGOMF[, c(2, 4)]))

showNum <- 6
pdf('MFbar.pdf', width = 15)
barplot(mf, showCategory = showNum)
dev.off()
pdf('MFdot.pdf', width = 15)
dotplot(mf, showCategory = showNum)
dev.off()

## convert to probe IDs
mf <- as.data.frame(mf)
mf[, 8] <- sapply(mf[, 8], function(x) {
  splitVec <- unlist(strsplit(x, split = '/', fixed = TRUE))
  convVec <- mergeProAnno[mergeProAnno[, 'KEGGID'] %in% splitVec, 'ProbeName']
  convVec <- paste(convVec, collapse = '/')
  return(convVec)
})

write.csv(mf, file = 'MFenrich.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~CC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
proGOCC <- proGO[proGO[, 3] == 'C', ]
cc <- enricher(DEGIDs, pvalueCutoff = 1, universe = unique(proGO[, 1]), minGSSize = 5, qvalueCutoff = 1, TERM2GENE = data.frame(proGOCC[, 2:1]), TERM2NAME = data.frame(proGOCC[, c(2, 4)]))

showNum <- 6
pdf('CCbar.pdf', width = 10)
barplot(cc, showCategory = showNum)
dev.off()
pdf('CCdot.pdf', width = 10)
dotplot(cc, showCategory = showNum)
dev.off()

## convert to probe IDs
cc <- as.data.frame(cc)
cc[, 8] <- sapply(cc[, 8], function(x) {
  splitVec <- unlist(strsplit(x, split = '/', fixed = TRUE))
  convVec <- mergeProAnno[mergeProAnno[, 'KEGGID'] %in% splitVec, 'ProbeName']
  convVec <- paste(convVec, collapse = '/')
  return(convVec)
})

write.csv(cc, file = 'CCenrich.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####################################################################
