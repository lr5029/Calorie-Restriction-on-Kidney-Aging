library(Seurat)
library(dplyr)
library(Matrix)
library(patchwork)

#read in the first sample
a1.data<-Read10X(data.dir ="./GSE137869_RAW/GSM4331828/") # this sub-folder contains GSM4331828_Kidney-M-Y_barcodes.tsv, GSM4331828_Kidney-M-Y_genes.tsv, and GSM4331828_Kidney-M-Y_matrix.mtx

# Error: raw_count shows "argument 'counts' is missing, with no default"
a1<-CreateSeuratObject(counts = a1.data, project = "M-Y") #this sample is Male, Young
#projects<-c("GSM4331829","GSM4331830","GSM4331831","GSM4331832","GSM4331833")
#sample_labels<-c("M-O","M-CR","F-Y","F-O","F-CR")

projects<-c("GSM4331829","GSM4331830","GSM4331831","GSM4331832","GSM4331833")
sample_labels<-c("M-O","M-CR","F-Y","F-O","F-CR")

#first read the raw data and create a list of objects
project_data<-vector("list", length = 5)
for (i in seq_along(projects)){
  temp.data<-Read10X(data.dir=paste("./GSE137869_RAW/",projects[i],sep="")) #read a new sample data
  project_data[i]<-CreateSeuratObject(counts = temp.data, project = sample_labels[i])
}

# Error: Could not find function "AddSamples"
# Error: Please provide a cell identifier for each object provided to merge
# for (i in seq_along(projects)){
#   project_data<-Read10X(data.dir=paste("./",projects[i],sep="")) #read a new sample data
#   a1 <- merge(x = a1, y = project_data, add.cell.ids = c("",sample_labels[i])) #add it to the existing data
# }

b1<-merge(x=a1,y=project_data,add.cell.ids=c("M-Y",sample_labels))

kidney<-b1
rm(a1)
rm(b1)

# percent of mitochondrial genes
kidney[["percent.mt"]] <- PercentageFeatureSet(kidney, pattern = "^Mt-")
VlnPlot(kidney, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot0 <- FeatureScatter(kidney, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot00 <- FeatureScatter(kidney, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot0 + plot00

# filter out unwanted cells
kidney <- subset(kidney, subset = nFeature_RNA > 200 & nFeature_RNA < 30000 & percent.mt < 20)

kidney <- NormalizeData(object = kidney, normalization.method = "LogNormalize", 
                    scale.factor = 10000)


pdf("./Graph/variable_genes.pdf")

# error: could not find function "FindVariableGenes"
# calculates the average expression and dispersion for each gene, places these genes into bins,
# and then calculates a z-score for dispersion within each bin.
kidney <- FindVariableFeatures(object = kidney, selection.method = "vst", nfeatures = 1050)

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(kidney), 20)

# plot variable features with labels
plot1 <- LabelPoints(plot = VariableFeaturePlot(kidney), points = top20, repel = TRUE)
print(plot1)
dev.off()

# Error: no slot of name "var.genes" for this object of class "Seurat"
# Error: no slot of name "RNA" for this object of class "Seurat"
length(x = kidney@assays[["RNA"]]@var.features)


# Scaling the data and removing unwanted sources of variation
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# Error: None of the requested variables to regress are present in the object.
kidney <- ScaleData(object = kidney, vars.to.regress = c("nCount_RNA", "percent.mt"))

# Perform linear dimensional reduction
kidney <- RunPCA(object = kidney, pc.genes = kidney@assays[["RNA"]]@var.features, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5)

pdf("./Graph/top_PC12_genes.pdf")
# Error: could not find function "VizPCA"
# VizPCA(object = kidney, pcs.use = 1:2)
plot2 <- VizDimLoadings(kidney, dims = 1:2, reduction = "pca")
print(plot2)
dev.off()

# visualizing both cells and genes that define the PCA
# Error: pdf has no content
pdf("./Graph/PCA.pdf")
# PCAPlot(object = kidney, dim.1 = 1, dim.2 = 2): error
plot3 <- DimPlot(kidney, reduction = "pca")
print(plot3)
dev.off()


# In particular DimHeatmap allows for easy exploration of the primary sources of
# heterogeneity in a dataset, and can be useful when trying to decide which PCs
# to include for further downstream analyses. Both cells and features are ordered
# according to their PCA scores
pdf("./Graph/PCA_heatmap.pdf")
# Error: PCHeatmap(object = kidney, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
plot4 <- DimHeatmap(kidney, dims = 1:6, cells = 500, balanced = TRUE)
print(plot4)
dev.off()


# Error: could not find function "PCElbowPlot"
# PCElbowPlot(object = kidney)
# PCA will create many components of variation but not all are useful. One way to
# select principal components that explain the variation in our data is to use an
# elbow or scree plot
ElbowPlot(kidney)

# Warning: The following arguments are not used: reduction.type, dims.use, print.output, save.SNN
# Suggested parameter: reduction instead of reduction.type; dims instead of dims.use; verbose instead of print.output

# First calculate k-nearest neighbors and construct the SNN graph (FindNeighbors), then run FindClusters
kidney <- FindNeighbors(kidney, reduction = "pca", dims = 1:10)

# The FindClusters function implements the procedure, and contains a resolution
# parameter that sets the 'granularity' of the downstream clustering, with increased
# values leading to a greater number of clusters.
kidney <- FindClusters(kidney, resolution = 0.2, verbose = 0, save.SNN = TRUE)

kidney <- RunTSNE(object = kidney, dims.use = 1:10, do.fast = TRUE, do.label= TRUE)


pdf("./Graph/tSNE_plot.pdf")
plot5 <- TSNEPlot(object = kidney, label=T)
print(plot5)
dev.off()

# Mapping gene expression onto cell clusters and assigning
# identity to clusters


# cluster 4
# pdf("./TSNE_Endo.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Nrp1", "Kdr"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # no match
# pdf("./TSNE_Podo.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Nphs1", "Nphs2"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # cluster 0 and 2 as one cluster
# pdf("./TSNE_PT.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Slc27a2", "Lrp2"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # cluster 10
# pdf("./TSNE_LOH.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Slc12a1", "Umod"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # little bit cluster 10
# pdf("./TSNE_DCT.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Slc12a3", "Pvalb"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # cluster 8
# pdf("./TSNE_CD-PC.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Aqp2", "Hsd11b2"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # cluster 5
# pdf("./TSNE_CD-1C.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Atp6v1g3", "Atp6v0d2"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # cluster 5
# pdf("./TSNE_CD-Trans.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Insrr", "Rhbg"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # no match
# pdf("./TSNE_Novel1.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Mki67", "Cdca3"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # cluster 6 and 7? 
# pdf("./TSNE_Fib.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Plac8", "S100a4"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # cluster 1
# pdf("./TSNE_Macro.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("C1qa", "C1qb"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # cluster 12
# pdf("./TSNE_Neutro.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("S100a8", "S100a9"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # no match
# pdf("./TSNE_B lymph.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Cd79a", "Cd79b"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # no match
# pdf("./TSNE_T lymph.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Ltb", "Cxcr6"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # cluster 6?
# pdf("./TSNE_NK.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = c("Gzma", "Nkg7"), cols = c("grey", "blue"), reduction = "tsne")
# dev.off()
# 
# 
# # no match
# pdf("./TSNE_Novel2.pdf", w=11, h=8.5)
# FeaturePlot(object = kidney, features = "Stmn1", cols = c("grey", "blue"), reduction = "tsne")
# dev.off()

# new.cluster.ids <- c("0", "1", "2", "3", "4", "5", "Endo",
#   "7", "8" , "9", "CD-PC", "11", "LOH", "13", "14", "15", "Neutro")

# names(new.cluster.ids) <- levels(kidney)
# kidney <- RenameIdents(kidney, new.cluster.ids)
# DimPlot(kidney, reduction = "tsne", label = TRUE, pt.size = 0.1) + NoLegend()

# pdf("./NewLabels.pdf", w=11, h=8.5)
# DimPlot(kidney, reduction = "tsne", label = TRUE, pt.size = 0.1) + NoLegend()
# dev.off()



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


# 0: MO = Macro; LEUK
# 1, 2: PT(S3); PT
# 4: Endo; Endo
# 11: Neutro; LEUK
# after looking at the link
# 3: LH(DL); diverse
# 5: IC-A; ICB or ICA
# 6: diverse; LEUK
# 7: diverse.
# 8: LH(AL); TAL
# 9: CD-PC; PC
# 10: MC; Fib or Mes
