########################KEGG ids#####################################
setwd('/home/Yulong/RESEARCH/BaiLinhan_array/Process/')

library('KEGGAPI')
library('NCBIAPI')


NCBIWholeMat <- read.table('ProteinTable62_22733.txt', sep = '\t', stringsAsFactors = FALSE)

mRNAidList <- geneid2mRNAid(NCBIWholeMat[, 6], n = 4)

mgrIDs <- getProID('mgr')
######################################################################
