#####Packages#####
if(!require(dplyr)){install.packages("dplyr")}
if(!require(Seurat)){install.packages("Seurat")}
if(!require(leiden)){install.packages("leiden")}
if(!require(clues)){install.packages("clues")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(ggpubr)){install.packages("ggpubr")}
if(!require(viridis)){install.packages("viridis")}
##########



#####Reference#####
#https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
##########



#####Input files#####
dir <- "/Users/kai/Desktop/convert_similarity_test/test_ICELL8"
folder <- "01_GBM"
convert <- "convert1.0.0"
vs <- "mappa0.9beta"

GBM.file.convert <- paste(dir, "/", folder, "/", convert, "_GBM.txt", sep="")
GBM.file.vs <- paste(dir, "/", folder, "/", vs, "_GBM.txt", sep="")
pdf.file <- paste(dir, "/", convert, "_vs_", vs, "_plots.pdf", sep="")
##########



#####Correlation#####
GBM.table.convert <-as.matrix(read.table(GBM.file.convert))
GBM.table.vs <- as.matrix(read.table(GBM.file.vs))
Cor.datapoints <- as.data.frame(cbind(as.vector(GBM.table.convert), as.vector(GBM.table.vs)))
r <- sprintf("%.2f", cor(Cor.datapoints[,1], Cor.datapoints[,2]))
Cor.datapoints <- Cor.datapoints %>% group_by_all() %>% summarise(COUNT = n())
Cor.datapoints$COUNT <- log10(Cor.datapoints$COUNT)
##########



#####Clustering#####
#process GBM.convert
GBM.seurat.convert <- CreateSeuratObject(read.table(GBM.file.convert))
GBM.seurat.convert <- NormalizeData(GBM.seurat.convert, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
GBM.seurat.convert <- FindVariableFeatures(GBM.seurat.convert, selection.method="vst", nfeatures=2000, verbose=FALSE)
GBM.seurat.convert <- ScaleData(GBM.seurat.convert, verbose=FALSE)
GBM.seurat.convert <- RunPCA(GBM.seurat.convert, features=VariableFeatures(object=GBM.seurat.convert))
GBM.seurat.convert <- RunUMAP(GBM.seurat.convert, reduction="pca", dims=1:50, umap.method="umap-learn", metric="correlation", seed.use=42)
GBM.seurat.convert <- FindNeighbors(GBM.seurat.convert, reduction="umap", dims=1:2)
GBM.seurat.convert <- FindClusters(GBM.seurat.convert, algorithm="leiden")
Cluster.convert <- as.data.frame(GBM.seurat.convert@reductions$umap@cell.embeddings)
Cluster.convert$Cluster <- GBM.seurat.convert@meta.data$seurat_clusters

#process GBM.vs
GBM.seurat.vs <- CreateSeuratObject(read.table(GBM.file.vs))
GBM.seurat.vs <- NormalizeData(GBM.seurat.vs, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
GBM.seurat.vs <- FindVariableFeatures(GBM.seurat.vs, selection.method="vst", nfeatures=2000, verbose=FALSE)
GBM.seurat.vs <- ScaleData(GBM.seurat.vs, verbose=FALSE)
GBM.seurat.vs <- RunPCA(GBM.seurat.vs, features=VariableFeatures(object=GBM.seurat.vs))
GBM.seurat.vs <- RunUMAP(GBM.seurat.vs, reduction="pca", dims=1:50, umap.method="umap-learn", metric="correlation")
GBM.seurat.vs <- FindNeighbors(GBM.seurat.vs, reduction="umap", dims=1:2)
GBM.seurat.vs <- FindClusters(GBM.seurat.vs, algorithm="leiden")
Cluster.vs <- as.data.frame(GBM.seurat.vs@reductions$umap@cell.embeddings)
Cluster.vs$Cluster <- GBM.seurat.vs@meta.data$seurat_clusters

#Adjusted Rand Index (ARI)
ARI <- sprintf("%.2f", as.numeric(adjustedRand(as.vector(Cluster.convert$Cluster), as.vector(Cluster.vs$Cluster))[1]))
##########



#####Generate plots#####



#plot corrilation
plot.Cor <- ggplot(Cor.datapoints, aes(x=V1, y=V2, color=COUNT)) +
    geom_abline(intercept = 0, slope = 1, color = "red", alpha=0.5) +
    geom_point() +
    theme(
        axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.85,0.25)
    ) +
    xlim(0, max(Cor.datapoints[,1:2])*1.05) +
    ylim(0, max(Cor.datapoints[,1:2])*1.05) +
    labs(
        x = paste("Count (", "convert1.0.0", ")", sep=""),
        y = paste("Count (", vs, ")", sep="")
    ) +
    scale_color_viridis(
        name = "Log10(count)"
    ) +
    geom_text(x = 0, y = max(Cor.datapoints[2]), hjust = 0, vjust = 1, label = paste("r=", r, sep=""), color="black")

#plot convert1.0.0 cluster
plot.cluster.convert <- ggplot(Cluster.convert, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) +
    geom_point() +
    theme(
        axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"
    ) +
    geom_text(x = min(Cluster.convert$UMAP_1), y = min(Cluster.convert$UMAP_2), hjust = 0, vjust = 0, label = "convert1.0.0", color = "black") +
    geom_text(x = max(Cluster.convert$UMAP_1), y = max(Cluster.convert$UMAP_2), hjust = 1, vjust = 1, label = paste("ARI=", ARI, sep=""), color = "black")

#plot vs cluster
plot.cluster.vs <- ggplot(Cluster.vs, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) +
    geom_point() +
    theme(
        axis.line = element_line(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"
    ) +
    geom_text(x = min(Cluster.vs$UMAP_1), y = min(Cluster.vs$UMAP_2), hjust = 0, vjust = 0, label = vs, color = "black")

#Combining plots
plot.total <-ggarrange(
    plot.Cor,
    plot.cluster.convert,
    plot.cluster.vs,
    ncol = 3, nrow = 1
    )
##########



#####Save results to .rda and .pdf file#####
#pdf file
pdf(file=pdf.file, width=12, height=4)
plot(plot.total)
dev.off()
##########
##########END##########
