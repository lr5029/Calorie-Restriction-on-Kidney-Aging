library(Seurat)
library(dplyr)
library(Matrix)
library(patchwork)
library(devtools)
# install.packages(c("WriteXLS", "devtools", "urltools", "pheatmap", "shinythemes","BiocManager"))
# install_github("yuliangwang/perturb.met")
# install.packages("remotes")
# remotes::install_github("ggjlab/scMCA")
# install_github("ggjlab/scMCA")
library(perturb.met)
library(scMCA)

cbPalette<-c("azure4","black","blue","brown","cadetblue","chartreuse","cyan",
             "darkorange","darkorchid","deeppink","gold","lightcoral","lightseagreen","magenta","red","lightsalmon","yellow","mediumorchid4","deepskyblue","mediumvioletred","olivedrab","cornsilk","lavender","navajowhite4")


load("kidney_aging_sc.RData")

kidney[["percent.mt"]] <- PercentageFeatureSet(kidney, pattern = "^Mt-")
VlnPlot(kidney, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#nFeature_RNA changed to 3000; percent.mt changed to 30. 11-10-2020
kidney <- subset(kidney, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 30)

kidney <- NormalizeData(object = kidney, normalization.method = "LogNormalize",
                        scale.factor = 10000)

pdf("./Graph/variable_genes.pdf")
kidney <- FindVariableFeatures(object = kidney, selection.method = "vst", nfeatures = 1000)
variable_genes_plot <- LabelPoints(plot = VariableFeaturePlot(kidney), points = top20, repel = TRUE)
print(variable_genes_plot)
dev.off()

kidney <- ScaleData(object = kidney, vars.to.regress = c("nCount_RNA", "percent.mt"))

# Perform linear dimensional reduction
kidney <- RunPCA(object = kidney, pc.genes = kidney@assays[["RNA"]]@var.features, do.print = TRUE, pcs.print = 1:5,
                 genes.print = 5)

pdf("./Graph/top_PC12_genes.pdf")
# Error: could not find function "VizPCA"
# VizPCA(object = kidney, pcs.use = 1:2)
top_12_genes_plot <- VizDimLoadings(kidney, dims = 1:2, reduction = "pca")
print(top_PC12_genes_plot)
dev.off()

pdf("./Graph/PCA.pdf")
# PCAPlot(object = kidney, dim.1 = 1, dim.2 = 2): error
pca_plot <- DimPlot(kidney, reduction = "pca")
print(pca_plot)
dev.off()

pdf("./Graph/PCA_heatmap.pdf")
# Error: PCHeatmap(object = kidney, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
heatmap_plot <- DimHeatmap(kidney, dims = 1:6, cells = 500, balanced = TRUE)
print(heatmap_plot)
dev.off()

ElbowPlot(kidney)

kidney <- FindNeighbors(kidney, reduction = "pca", dims = 1:10)

kidney <- FindClusters(kidney, resolution = 0.2, verbose = 0, save.SNN = TRUE)


kidney <- RunUMAP(kidney, dims = 1:10)

pdf("./Graph/dim_plot.pdf")
plot <- DimPlot(kidney, reduction = "umap",cols=cbPalette[1:16], label = TRUE)
print(plot)
dev.off()


MIN_LOGFOLD_CHANGE = 2 # set to minimum required average log fold change in gene expression.
MIN_PCT_CELLS_EXPR_GENE = .25  # minimum percent of cells that must express gene in either clstr.

# Here we find all markers for 
all.markers = FindAllMarkers(kidney,
                             min.pct = MIN_PCT_CELLS_EXPR_GENE,
                             logfc.threshold = MIN_LOGFOLD_CHANGE,
                             only.pos = TRUE)
marker_genes_cluster <- all.markers %>%
  group_by(cluster) %>%
  top_n(20, avg_logFC)

# Based on http://humphreyslab.com/SingleCell/search.php and scMCA
# 0: PT(S1-S2)
# 1: Macrophage_Ccl4 high (Kidney)
# 2: PT(S3)
# 3: IC-A
# 4: LH(AL)
# 5: LH(DL)
# 6: EC
# 7: T cell(Kidney)
# 8: Macrophage_Lyz2 high(Kidney)
# 9: CD-PC
# 10: MC
# 11: Neutrophil progenitor_S100a8 high(Kidney)

# Use scMCA to identify cluster 7,8,11 since using Kit database cannot identify
# a identical cluster

cluster_data<-as.matrix(kidney@assays$RNA@counts)
total<-colSums(cluster_data)
total<-total/median(total)
memory.limit(size=56000)
cluster_data<-cluster_data/matrix(rep(total,nrow(cluster_data)),nrow=nrow(cluster_data),ncol=ncol(cluster_data),byrow = T)
top1000<-head(VariableFeatures(kidney), 1000)
cluster_data<-cluster_data[rownames(cluster_data) %in% top1000,]
mca_result<-scMCA(scdata = cluster_data,numbers_plot = 3)
head(rev(sort(table(mca_result$scMCA[kidney$seurat_clusters==0]))/sum(kidney$seurat_clusters==0)))
head(rev(sort(table(mca_result$scMCA[kidney$seurat_clusters==8]))/sum(kidney$seurat_clusters==8)))
head(rev(sort(table(mca_result$scMCA[kidney$seurat_clusters==11]))/sum(kidney$seurat_clusters==11)))


new.cluster.ids <- 
  c("PT(S1-S2)", "Macrophage_Ccl4 high", "PT(S3)", "IC-A", "LH(AL)",
    "LH(DL)", "EC", "T cell", "Macrophage_Lyz2 high", "CD-PC", "MC", "Neutrophil progenitor_S100a8 high")

names(new.cluster.ids) <- levels(kidney)
kidney <- RenameIdents(kidney, new.cluster.ids)
pdf("./Graph/dim_plot_label.pdf")
plot_label <- DimPlot(kidney, reduction = "umap", label = TRUE, pt.size = 0.1)
print(plot_label)
dev.off()
