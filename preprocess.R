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

save(normAveArray, file = '/home/Yulong/RESEARCH/BaiLinhan_array/Process/normAveArray.RData')
######################################################################


################################compare with GEO platform##############
library('limma')

load('/home/Yulong/RESEARCH/BaiLinhan_array/Process/normAveArray.RData')

GPL9023Anno <- read.csv('/home/Yulong/RESEARCH/BaiLinhan_array/Process/GPL9023_anno.txt', sep = '\t', comment.char = '#', stringsAsFactor = FALSE)

rawAnno <- normAveArray$genes

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~inter probes~~~~~~~~~~~~~~~~~~~~
rawAnnoProbes <- rawAnno[, 'ProbeName']
GPL9023AnnoProbes <- GPL9023Anno[, 'ID']

interProbes <- intersect(rawAnnoProbes, GPL9023AnnoProbes)

combineProbes <- cbind(GPL9023Anno[match(interProbes, GPL9023AnnoProbes), ],
              rawAnno[match(interProbes, rawAnnoProbes), ])

sum(combineProbes[, 'SEQUENCE'] == combineProbes[, 'Sequence'])
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~select Magnaporthe grisea~~~~~~~~~~~~~~~~~~~
mgrP <- GPL9023Anno[GPL9023Anno[, 'ORGANISM'] == 'Magnaporthe grisea', 'ID']
mgrPIdx <- match(mgrP, rownames(normAveArray$E))

mgrArray <- cbind(normAveArray$E[mgrPIdx, ],
                  normAveArray$genes[mgrPIdx, ])

write.csv(mgrArray, file = '/home/Yulong/RESEARCH/BaiLinhan_array/Process/mgrArray.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#######################################################################
