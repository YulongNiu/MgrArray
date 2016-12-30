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

##~~~~~~~~~~~~~~~~~~~~KEGG ID to NCBI gene ID~~~~~~~~~~~~~~~~~~~~~
KEGGIDMat <- convKEGG('ncbi-geneid', 'mgr')
KEGGIDMat[, 2] <- sapply(strsplit(KEGGIDMat[, 2], split = ':', fixed = TRUE), '[[', 2)
colnames(KEGGIDMat) <- c('KEGGID', 'geneid')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~merge IDs with array probes~~~~~~~~~~~~~~~~~~
load('/home/Yulong/RESEARCH/BaiLinhan_array/Process/mgrArray.RData')

mergeAnno <- merge(KEGGIDMat, mRNAMat, by.x = 'geneid', by.y = 'geneid')

## with array probes
arrayAnno <- mgrArray$genes[, c('ProbeName', 'RefSeqAccession')]
mergeGeneAnno <- merge(mergeAnno, arrayAnno, by.x = 'mRNAid', by.y = 'RefSeqAccession')
save(mergeGeneAnno, file = 'mergeGeneAnno.RData')


## with ncbi-proteinid and uniprot
proidMat <- convKEGG('ncbi-proteinid', 'mgr')
proidMat[, 2] <- sapply(strsplit(proidMat[, 2], split = ':', fixed = TRUE), '[[', 2)
colnames(proidMat) <- c('KEGGID', 'proteinid')

unidMat <- convKEGG('uniprot', 'mgr')
unidMat[, 2] <- sapply(strsplit(unidMat[, 2], split = ':', fixed = TRUE), '[[', 2)
colnames(unidMat) <- c('KEGGID', 'uniprotid')

mergeProUniMat <- merge(proidMat, unidMat, by.x = 'KEGGID', by.y = 'KEGGID')
mergeProAnno <- merge(mergeGeneAnno, mergeProUniMat, by.x = 'KEGGID', by.y = 'KEGGID')
save(mergeProAnno, file = 'mergeProAnno.RData')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################################


