#######################map ID###########################
## map ID from http://pevsnerlab.kennedykrieger.org/flatfiles/Mgr.data-output-5.txt
standardMat <- read.table('/home/Yulong/RESEARCH/BaiLinhan_array/Process/Mgr.data-output-5.txt', stringsAsFactors = FALSE)

mgrArray <- read.csv('/home/Yulong/RESEARCH/BaiLinhan_array/Process/mgrArray.csv', row.names = 1, stringsAsFactors = FALSE)

hasLogic <- mgrArray[, 'GeneName'] %in% standardMat[, 2]
head(mgrArray[!hasLogic, ])
########################################################
