##########################limma DEGs############################
library('limma')
library('pheatmap')
library('ggplot2')
library('directlabels')

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
topNum <- 100
topLogic <- rownames(yeastArray) %in% rownames(DEGMat)[1:topNum]
heatmapCount <- yeastArray$E[topLogic, ]

annoCol <- data.frame(Group = rep(c('Blank', 'Yeast'), each = 3))
row.names(annoCol) <- colnames(heatmapCount)
annoColor <- list(Group = c(Blank = '#00C19F', Yeast = '#F8766D'))

pdf('/home/Yulong/RESEARCH/BaiLinhan_array/Process/heatmap_100.pdf')
pheatmap(heatmapCount, annotation_col = annoCol, annotation_colors = annoColor, fontsize = 8, fontsize_row = 4, annotation_legend = TRUE)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~
pca <- prcomp(t(yeastArray$E))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = rep(c('Blank', 'Yeast'), each = 3), ID = colnames(yeastArray))

pdf('/home/Yulong/RESEARCH/BaiLinhan_array/Process/PCA.pdf')
groupCol <- c('#00C19F', '#F8766D')
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = groupCol) +
  geom_dl(aes(label = ID, color = Group), method = 'smart.grid')
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################
