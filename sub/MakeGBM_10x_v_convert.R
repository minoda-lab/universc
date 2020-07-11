#####Packages#####
if(!require(Seurat)){install.packages("Seurat")}
##########



#####Input files#####
Cellranger.INPUT <- "tiny-count-v3-v302/outs/outs/filtered_feature_bc_matrix"
Convert.INPUT <- "test-10x-v3-v302/outs/filtered_feature_bc_matrix"
Out.dir <- "test_10x_GBM"
##########

#####Read in input files#####
Cellranger.GBM <- as.data.frame(as.matrix(Read10X(data.dir=Cellranger.INPUT, gene.column = 1)))
Convert.GBM <- as.data.frame(as.matrix(Read10X(data.dir=Convert.INPUT, gene.column = 1)))
##########

#####Sort and clean GBM#####
Cellranger.GBM <- Cellranger.GBM[order(rownames(Cellranger.GBM)),order(colnames(Cellranger.GBM))]
for (i in 1:ncol(Cellranger.GBM)) {
    colnames(Cellranger.GBM)[i] <- strsplit(colnames(Cellranger.GBM)[i], "-", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1]
}
Convert.GBM <- Convert.GBM[order(rownames(Convert.GBM)),order(colnames(Convert.GBM))]
for (i in 1:ncol(Convert.GBM)) {
    colnames(Convert.GBM)[i] <- strsplit(colnames(Convert.GBM)[i], "-", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1]
}
##########

#####Matching barcodes#####
barcode.intercept <- colnames(Cellranger.GBM)[which(colnames(Cellranger.GBM)%in%colnames(Convert.GBM))]
Cellranger.GBM <- Cellranger.GBM[,which(colnames(Cellranger.GBM)%in%barcode.intercept)]
Convert.GBM <- Convert.GBM[,which(colnames(Convert.GBM)%in%barcode.intercept)]
##########

#####Matching genes#####
gene.union <- unique(c(rownames(Cellranger.GBM), rownames(Convert.GBM)))[order(unique(c(rownames(Cellranger.GBM), rownames(Convert.GBM))))]
for (i in 1:length(gene.union)) {
    if (!gene.union[i]%in%rownames(Cellranger.GBM)) {
        Cellranger.GBM[nrow(Cellranger.GBM) + 1,] <- rep(0, ncol(Cellranger.GBM))
        rownames(Cellranger.GBM)[nrow(Cellranger.GBM)] <- gene.union[i]
    }
}
Cellranger.GBM <- Cellranger.GBM[order(rownames(Cellranger.GBM)),]
for (i in 1:length(gene.union)) {
    if (!gene.union[i]%in%rownames(Convert.GBM)) {
        Convert.GBM[nrow(Convert.GBM) + 1,] <- rep(0, ncol(Convert.GBM))
        rownames(Convert.GBM)[nrow(Convert.GBM)] <- gene.union[i]
    }
}
Convert.GBM <- Convert.GBM[order(rownames(Convert.GBM)),]
##########

#####Output GBM#####
system(paste("mkdir -p",  Out.dir)
Cellranger.file <- paste(Out.dir, "cellranger3.0.2_GBM.txt", sep="/")
Convert.file <- paste(Out.dir, "convert1.0.0_GBM.txt", sep="/")

write.table(Cellranger.GBM, file=Cellranger.file, quote=FALSE,sep="\t", eol="\n", row.names=TRUE, col.names=NA)
write.table(Convert.GBM, file=Convert.file, quote=FALSE,sep="\t", eol="\n", row.names=TRUE, col.names=NA)
##########

##########END##########
