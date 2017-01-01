########################convert GO##################################
library('clusterProfiler')
library('stringr')
library('GO.db')

load('mergeProAnno.RData')
mergeProAnno[, 1] <- sapply(strsplit(as.character(mergeProAnno[, 1]), split = ':', fixed = TRUE), '[[', 2)
DEGMat <- read.csv('DEGMat.csv', stringsAsFactor = FALSE, row.names = 1)
load('mergeGeneAnno.RData')

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

## BP
proGOBP <- proGO[proGO[, 3] == 'P', ]
bp <- enricher(DEGIDs, pvalueCutoff = 1, universe = unique(proGO[, 1]), minGSSize = 5, qvalueCutoff = 1, TERM2GENE = data.frame(proGOBP[, 2:1]), TERM2NAME = data.frame(proGOBP[, c(2, 4)]))
write.csv(bp, file = 'BPenrich.csv')

####################################################################
