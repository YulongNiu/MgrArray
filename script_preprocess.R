###############################preprocess array#######################
library('limma')

## read in raw array
storePath <- '/home/Yulong/RESEARCH/BaiLinhan_array/Array/oridata'
filePaths <- dir(storePath)
## targets <- data.frame(SampleNumber = 1:length(filePaths),
##                       FileNamae = filePaths,
##                       Conditions = paste0('array', 1:length(filePaths)))
arrayNames <- c(paste0('YeastExtract', 1:3), paste0('Blank', 1:3), paste0('EucalyptusOil', 1:3), paste0('Chemical', 1:3))
rawArray <- read.maimages(filePaths, source = 'agilent', green.only = TRUE, path = storePath, names = arrayNames)


## background correct
bgArray <- backgroundCorrect(rawArray, method = 'normexp', offset=50)

## normalization
normArray <- normalizeBetweenArrays(bgArray, method = 'quantile')
normAveArray <- avereps(normArray, ID=normArray$genes$ProbeName)

save(normAveArray, file = '/home/Yulong/RESEARCH/BaiLinhan_array/Process/normAveArray.RData')
######################################################################


################################compare with GEO platform##############
library('limma')

load('/home/Yulong/RESEARCH/BaiLinhan_array/Process/normAveArray.RData')

GPL9203Anno <- read.csv('/home/Yulong/RESEARCH/BaiLinhan_array/Process/GPL9203_anno.txt', sep = '\t', comment.char = '#', stringsAsFactor = FALSE)

singleRaw <- read.csv('/home/Yulong/RESEARCH/BaiLinhan_array/Process/single_array_raw.csv', stringsAsFactor = FALSE)

rawAnno <- normAveArray$genes

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~inter probes~~~~~~~~~~~~~~~~~~~~
rawAnnoProbes <- rawAnno[, 'ProbeName']
GPL9203AnnoProbes <- GPL9203Anno[, 'ID']

interProbes <- intersect(rawAnnoProbes, GPL9203AnnoProbes)

combineProbes <- cbind(GPL9203Anno[match(interProbes, GPL9203AnnoProbes), ],
              rawAnno[match(interProbes, rawAnnoProbes), ])

sum(combineProbes[, 'SEQUENCE'] == combineProbes[, 'Sequence'])
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~select Magnaporthe grisea~~~~~~~~~~~~~~~~~~~
mgrP <- GPL9203Anno[GPL9203Anno[, 'ORGANISM'] == 'Magnaporthe grisea', 'ID']
mgrPIdx <- match(mgrP, rownames(normAveArray$E))

mgrArray <- normAveArray[mgrPIdx, ]

## add anotation
mgrAddAnno <- singleRaw[match(mgrP, singleRaw[, 1]), 38:54]
mgrArray$genes <- cbind(mgrArray$genes, mgrAddAnno)

save(mgrArray, file = '/home/Yulong/RESEARCH/BaiLinhan_array/Process/mgrArray.RData')

mgrArrayData <- cbind(mgrArray$E,
                      mgrArray$genes)

write.csv(mgrArrayData, file = '/home/Yulong/RESEARCH/BaiLinhan_array/Process/mgrArrayData.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################


########################KEGG ids######################################
library('KEGGAPI')

mgrIDs <- getProID('mgr')
######################################################################
