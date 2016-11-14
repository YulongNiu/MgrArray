##########################limma DEGs############################
library('limma')

load('/home/Yulong/RESEARCH/BaiLinhan_array/Process/mgrArray.RData')

yeastArray <- mgrArray[, 6:1]
group <- rep(c('blank', 'yeast'), each = 3)
design <- model.matrix(~factor(group))
colnames(design) <- c('blank', 'yeast+vs')

fit <- lmFit(yeastArray, design)
fit <- eBayes(fit)
DEGMat <- topTable(fit, coef = 2, n = nrow(yeastArray), adjust = "BH")
write.csv(DEGMat[, c(-29, -30, -33)], file = '/home/Yulong/RESEARCH/BaiLinhan_array/Process/DEGMat.csv')

##~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################
