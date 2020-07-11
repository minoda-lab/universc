#####Packages#####
if(!require(Seurat)){install.packages("Seurat")}
##########



#####Input files#####
Nadia.INPUT <- "/Users/kai/Desktop/convert_similarity_test/test_Nadia/INPUT/dropseqpipe0.6/samples/SRR1873277_S1_L001/umi"
Convert.INPUT <- "/Users/kai/Desktop/convert_similarity_test/test_Nadia/INPUT/convert1.0.0/outs/filtered_feature_bc_matrix"
##########

#####Read in input files#####
Nadia.GBM <- as.data.frame(as.matrix(Read10X(data.dir=Nadia.INPUT, gene.column = 1)))
Convert.GBM <- as.data.frame(as.matrix(Read10X(data.dir=Convert.INPUT, gene.column = 2)))
##########

#####Sort and clean GBM#####
Nadia.GBM <- Nadia.GBM[order(rownames(Nadia.GBM)),order(colnames(Nadia.GBM))]
Convert.GBM <- Convert.GBM[order(rownames(Convert.GBM)),order(colnames(Convert.GBM))]
for (i in 1:ncol(Convert.GBM)) {
    colnames(Convert.GBM)[i] <- strsplit(colnames(Convert.GBM)[i], "-", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1]
}
##########

#####Matching barcodes#####
barcode.intercept <- colnames(Nadia.GBM)[which(colnames(Nadia.GBM)%in%colnames(Convert.GBM))]
Nadia.GBM <- Nadia.GBM[,which(colnames(Nadia.GBM)%in%barcode.intercept)]
Convert.GBM <- Convert.GBM[,which(colnames(Convert.GBM)%in%barcode.intercept)]
##########

#####Matching genes#####
gene.union <- unique(c(rownames(Nadia.GBM), rownames(Convert.GBM)))[order(unique(c(rownames(Nadia.GBM), rownames(Convert.GBM))))]
for (i in 1:length(gene.union)) {
    if (!gene.union[i]%in%rownames(Nadia.GBM)) {
        Nadia.GBM[nrow(Nadia.GBM) + 1,] <- rep(0, ncol(Nadia.GBM))
        rownames(Nadia.GBM)[nrow(Nadia.GBM)] <- gene.union[i]
    }
}
Nadia.GBM <- Nadia.GBM[order(rownames(Nadia.GBM)),]
for (i in 1:length(gene.union)) {
    if (!gene.union[i]%in%rownames(Convert.GBM)) {
        Convert.GBM[nrow(Convert.GBM) + 1,] <- rep(0, ncol(Convert.GBM))
        rownames(Convert.GBM)[nrow(Convert.GBM)] <- gene.union[i]
    }
}
Convert.GBM <- Convert.GBM[order(rownames(Convert.GBM)),]
##########

#####Output GBM#####
Out.dir <- "/Users/kai/Desktop/convert_test/test_Nadia/01_GBM"
Nadia.file <- paste(Out.dir, "dropseqpipe0.6_GBM.txt", sep="/")
Convert.file <- paste(Out.dir, "convert1.0.0_GBM.txt", sep="/")

write.table(Nadia.GBM, file=Nadia.file, quote=FALSE,sep="\t", eol="\n", row.names=TRUE, col.names=NA)
write.table(Convert.GBM, file=Convert.file, quote=FALSE,sep="\t", eol="\n", row.names=TRUE, col.names=NA)
##########

##########END##########
