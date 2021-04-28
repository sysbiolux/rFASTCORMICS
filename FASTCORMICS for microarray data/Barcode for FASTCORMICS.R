#-----------------------#
# Set working directory #
#-----------------------#

# First, set the path to the folder containing the .CEL files

setwd("PATH TO .CEL FILES") # FASTCORMICS for microarray data/ExampleData/CEL example files (from GSE49910)


#------------------#
# Install packages #
#------------------#

install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hugene.1.0.st.v1frmavecs")
BiocManager::install("hgu133afrmavecs")
BiocManager::install("hgu133plus2frmavecs")
BiocManager::install("hgu133a2frmavecs")

BiocManager::install("oligo")

BiocManager::install("frma")

BiocManager::install("affy")


#--------------#
# load vectors #
#--------------#
# available vectors can be found at http://barcode.luhs.org/index.php?page=intro
# vectors containing median and standard deviation of the mode corresponding to the unexpressed genes
# choose which one corresponds to your data

## for hugene.1.0.st.v1
library("hugene.1.0.st.v1frmavecs")
data("hugene.1.0.st.v1frmavecs")
data("hugene.1.0.st.v1barcodevecs")

## for hgu133a
library("hgu133afrmavecs")
data("hgu133afrmavecs")
data("hgu133abarcodevecs")

## for hgu133plus2
library("hgu133plus2frmavecs")
data("hgu133plus2frmavecs")
data("hgu133plus2barcodevecs")

## for hgu133a2
library("hgu133a2frmavecs")
data("hgu133a2frmavecs")
data("hgu133a2barcodevecs")

#----------------#
# read cel files #
#----------------#

# library("oligo")
# celFiles = list.celfiles() #in current folder
# rawFiles = read.celfiles(celFiles)

library("affy")
celFiles = list.celfiles() #in current folder
rawFiles = read.affybatch(celFiles)


#-------------------------#
# Normalization & BARCODE #
#-------------------------#

# frma normalize the data and retrieve z-scores via Barcode
# The Barcode 3.0 paper has been published. Please view it here: McCall et al., 
# "The Gene Expression Barcode 3.0: improved data processing and mining tools"
# Nucleic Acids Research. 2014 Jan;42:D938-D943.  

# make sure to use the correct vectors!

library("frma")
rawFiles.eset <- frma(rawFiles, summarize="median_polish", target="core", input.vecs=hgu133plus2frmavecs)

rawFiles.barcode <- barcode(exprs(rawFiles.eset), platform=NULL, mu=hgu133plus2barcodevecs$mu, 
                            tau=hgu133plus2barcodevecs$tau, output="z-score")

#-----------#
# save data #
#-----------#

write.table(rawFiles.barcode, "barcode.txt", quote=F, col.names=F, row.names=F)
write.table(colnames(rawFiles.barcode),"colnames.txt", quote=F, col.names=F, row.names=F)
write.table(rownames(rawFiles.barcode),"rownames.txt", quote=F, col.names=F, row.names=F)
write.table(exprs(rawFiles.eset), "frma.txt", quote=F, col.names=F, row.names=F)
