setwd("P:/PhD Project/2019_Jessica_Ysaline/GSE46824 dataset/GSE46824/CEL unzipped")

library("oligo")
library(pd.hugene.1.0.st.v1)

celFiles = list.celfiles()
rawFiles = read.celfiles(celFiles)

library("hugene.1.0.st.v1frmavecs")# vector containing median and standard deviation of the mode corresponding to the unexpressed genes
data("hugene.1.0.st.v1frmavecs")
data("hugene.1.0.st.v1barcodevecs")


library("frma")
rawFiles.eset <- frma(rawFiles, summarize="median_polish", target="core", input.vecs=hugene.1.0.st.v1frmavecs)

rawFiles.barcode <- barcode(exprs(rawFiles.eset), platform=NULL, mu=hugene.1.0.st.v1barcodevecs$mu, 
                            tau=hugene.1.0.st.v1barcodevecs$tau, output="z-score")

# write.table(rawFiles.barcode, "barcode_GSE46824.txt", quote=F, col.names=F, row.names=F)
# write.table(colnames(rawFiles.barcode),"col_GSE46824.txt", quote=F, col.names=F, row.names=F)
# write.table(rownames(rawFiles.barcode),"row_GSE46824.txt", quote=F, col.names=F, row.names=F)
# write.table(exprs(rawFiles.eset), "frma_GSE46824.txt", quote=F, col.names=F, row.names=F)

library("ff")
rawFiles.rma <- as.data.frame(as.ffdf(exprs(rawFiles[,])))
# write.table(rawFiles.rma, "rma_GSE46824.txt", quote=F, col.names=F, row.names=F)


# Quality control
library("affy")
library("gplots")

boxplot(log2(rawFiles.rma),main="Before fRMA normalization (GSE46824)", xlab="Microarray names", ylab="Intensity (log2)")
plotDensity(log2(rawFiles.rma))
rawFiles.rma.cor = cor(log2(rawFiles.rma), method="pearson")
a=heatmap.2(rawFiles.rma.cor, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.6, cexCol=0.6, margins=c(11,11))


boxplot(exprs(rawFiles.eset),main="After fRMA normalization (GSE46824)", xlab="Microarray names", ylab="Intensity (log2)")
plotDensity(exprs(rawFiles.eset))
rawFiles.eset.cor = cor(exprs(rawFiles.eset), method="pearson")
a=heatmap.2(rawFiles.eset.cor, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.6, cexCol=0.6, margins=c(11,11))

boxplot(rawFiles.barcode,main="Barcode (GSE46824)", xlab="Microarray names", ylab="z-scores")
plotDensity(rawFiles.barcode)
rawFiles.barcode = cor(rawFiles.barcode, method="pearson")
a=heatmap.2(rawFiles.barcode, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.6, cexCol=0.6, margins=c(11,11))


