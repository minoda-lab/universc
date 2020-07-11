#####Packages#####
if(!require(Seurat)){install.packages("Seurat")}
##########



#####Input files#####
ICELL8.INPUT <- "/Users/kai/Desktop/convert_similarity_test/test_ICELL8/INPUT/mappa0.9beta/test-icell8-custom-mappa-0.9_genematrix.csv"
Convert.INPUT <- "/Users/kai/Desktop/convert_similarity_test/test_ICELL8/INPUT/convert1.0.0/outs/filtered_feature_bc_matrix"
##########

#####Read in input files#####
ICELL8.GBM <- read.csv(ICELL8.INPUT)
Convert.GBM <- as.data.frame(as.matrix(Read10X(data.dir=Convert.INPUT, gene.column = 1)))
##########

#####Sort and clean GBM#####
rownames(ICELL8.GBM) <- ICELL8.GBM[,1]
ICELL8.GBM <- ICELL8.GBM[,-1]
ICELL8.GBM <- ICELL8.GBM[order(rownames(ICELL8.GBM)),order(colnames(ICELL8.GBM))]
Convert.GBM <- Convert.GBM[order(rownames(Convert.GBM)),order(colnames(Convert.GBM))]
for (i in 1:ncol(Convert.GBM)) {
    colnames(Convert.GBM)[i] <- strsplit(colnames(Convert.GBM)[i], "-", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1]
}
##########

#####Matching barcodes#####
barcode.intercept <- colnames(ICELL8.GBM)[which(colnames(ICELL8.GBM)%in%colnames(Convert.GBM))]
ICELL8.GBM <- ICELL8.GBM[,which(colnames(ICELL8.GBM)%in%barcode.intercept)]
Convert.GBM <- Convert.GBM[,which(colnames(Convert.GBM)%in%barcode.intercept)]
##########

#####Matching genes#####
gene.union <- unique(c(rownames(ICELL8.GBM), rownames(Convert.GBM)))[order(unique(c(rownames(ICELL8.GBM), rownames(Convert.GBM))))]
for (i in 1:length(gene.union)) {
    if (!gene.union[i]%in%rownames(ICELL8.GBM)) {
        ICELL8.GBM[nrow(ICELL8.GBM) + 1,] <- rep(0, ncol(ICELL8.GBM))
        rownames(ICELL8.GBM)[nrow(ICELL8.GBM)] <- gene.union[i]
    }
}
ICELL8.GBM <- ICELL8.GBM[order(rownames(ICELL8.GBM)),]
for (i in 1:length(gene.union)) {
    if (!gene.union[i]%in%rownames(Convert.GBM)) {
        Convert.GBM[nrow(Convert.GBM) + 1,] <- rep(0, ncol(Convert.GBM))
        rownames(Convert.GBM)[nrow(Convert.GBM)] <- gene.union[i]
    }
}
Convert.GBM <- Convert.GBM[order(rownames(Convert.GBM)),]
##########

#####Output GBM#####
Out.dir <- "/Users/kai/Desktop/convert_similarity_test/test_ICELL8/01_GBM"
ICELL8.file <- paste(Out.dir, "mappa0.9beta_GBM.txt", sep="/")
Convert.file <- paste(Out.dir, "convert1.0.0_GBM.txt", sep="/")

write.table(ICELL8.GBM, file=ICELL8.file, quote=FALSE,sep="\t", eol="\n", row.names=TRUE, col.names=NA)
write.table(Convert.GBM, file=Convert.file, quote=FALSE,sep="\t", eol="\n", row.names=TRUE, col.names=NA)
##########

##########END##########
