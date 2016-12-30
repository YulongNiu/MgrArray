########################KEGG ids#####################################
##~~~~~~~~~~~~~~~~~~~~~~~NCBI gene id to mRNA accession num~~~~~~~~~~~~
setwd('/home/Yulong/RESEARCH/BaiLinhan_array/Process/')

library('KEGGAPI')
library('NCBIAPI')
library('ParaMisc')


NCBIWholeMat <- read.csv('ProteinTable62_22733.txt', sep = '\t', stringsAsFactors = FALSE)

cutMat <- CutSeqEqu(nrow(NCBIWholeMat), 200)
saveFolder <- '/home/Yulong/RESEARCH/BaiLinhan_array/Process/saveNCBIid/'

for (i in 1:ncol(cutMat)) {
  mRNAidList <- geneid2mRNAid(NCBIWholeMat[cutMat[1, i] : cutMat[2, i], 6], n = 4)
  fileName <- paste0(saveFolder, 'mRNAid_', cutMat[1, i], '_', cutMat[2, i], '.RData')
  save(mRNAidList, file = fileName)
}

NCBIfiles <- dir('saveNCBIid', full.names = TRUE)
mRNAMat <- matrix(ncol = 2)

for (i in NCBIfiles) {
  load(i)
  eachLen <- sapply(mRNAidList, length)
  eachMat <- cbind(rep(names(mRNAidList), eachLen),
                   unlist(mRNAidList))
  mRNAMat <- rbind(mRNAMat, eachMat)
}

mRNAMat <- mRNAMat[-1, ]
rownames(mRNAMat) <- NULL
colnames(mRNAMat) <- c('geneid', 'mRNAid')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~KEGG ID to NCBI gene ID~~~~~~~~~~~~~~~~~~~~~~~~
KEGGIDMat <- convKEGG('ncbi-geneid', 'mgr')
KEGGIDMat[, 2] <- sapply(strsplit(KEGGIDMat[, 2], split = ':', fixed = TRUE), '[[', 2)
colnames(KEGGIDMat) <- c('KEGGID', 'geneid')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## merge IDs
mergeAnno <- merge(KEGGIDMat, mRNAMat, by.x = 'geneid', by.y = 'geneid')
save(mergeAnno, file = 'mergeAnno.RData')
######################################################################
