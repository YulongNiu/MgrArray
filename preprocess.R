###############################preprocess array#######################
library('limma')

## read in raw array
storePath <- '/home/Yulong/RESEARCH/BaiLinhan_array/Array/oridata'
filePaths <- dir(storePath)
## targets <- data.frame(SampleNumber = 1:length(filePaths),
##                       FileNamae = filePaths,
##                       Conditions = paste0('array', 1:length(filePaths)))
rawArray <- read.maimages(filePaths, source = 'agilent', green.only = TRUE, path = storePath, names = paste0('array', 1:length(filePaths)))


## background correct
bgArray <- backgroundCorrect(rawArray, method = 'normexp', offset=50)

## normalization
normArray <- normalizeBetweenArrays(bgArray, method = 'quantile')
normAveArray <- avereps(normArray, ID=normArray$genes$ProbeName)

######################################################################
