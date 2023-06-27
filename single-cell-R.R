#loading in the libraries that are required for scRNA-seq analysis 
#using Seurat for (QC, normalization, scaling, features selection, cluster)
#using scType for cell annotation

library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(DropletUtils)
library(patchwork)
library(harmony)
library(SoupX)
library(glmnet)
library(ggplot2)
library(ComplexHeatmap)

#for cell annotation
library(SCINA)
library(harmony)
library(cerebroApp)
library(msigdbi)
library(scmap)
library(celldex)
library(SingleR)


##### Below is the standard flow that is executed #https://satijalab.org/seurat/articles/essential_commands.html
pbmc.counts <- Read10X(data.dir = "counts/")
pbmc <- CreateSeuratObject(counts = pbmc.counts)
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc)
pbmc <- RunTSNE(object = pbmc)
DimPlot(object = pbmc, reduction = "tsne")


##scRNA-seq analysis pipeline

#Not just 10x data, you can read any scRNA-seq technology quantification data 

#raw_data <- read.table("/path/to/file.tsv")
#object <- CreateSeuratObject(counts = raw_data, ...)


#################For a single phenotype data##################

#to read 10x data specifically
tumor_data = "tumor.out"
tumor_data = Read10X(data.dir = tumor_data)  # Seurat function to read in 10x count data
dim(tumor_data)

pc = CreateSeuratObject(tumor_data, project = "Prostate Tumor")

pc@meta.data$disease = "tumor_prostate_gland"
pc@meta.data$orig.ident.

#Quality Control (identification of dead cells and subsequent removal)

vln_plot1 <- VlnPlot(pc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
vln_plot1
ggsave("feature_count_RNA.png", vln_plot1, dpi = 600)

pc[["percent.mt"]] <- PercentageFeatureSet(pc, pattern = "^MT-")
vln_plot2 <- VlnPlot(pc, features = c("percent.mt"), ncol = 1)
ggsave("percent_mt.png", vln_plot2, dpi = 600)

#filtering dead cells out, highly expressive cells and lowly expressive cells are being filtered out
##nFeature_RNA is the number of genes detected in each cell. 
##nCount_RNA is the total number of molecules detected within a cell. 
##Low nFeature_RNA for a cell indicates that it may be dead/dying or an empty droplet. 
##High nCount_RNA and/or nFeature_RNA indicates that the "cell" may in fact be a doublet (or multiplet). 
##In combination with %mitochondrial reads, removing outliers from these groups removes most doublets/dead cells/empty droplets, hence why filtering is a common pre-processing step.
pc = subset(pc, subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.mt < 10)

pc[["percent.mt"]] <- PercentageFeatureSet(pc, pattern = "^MT-")
vln_plot3 <- VlnPlot(pc, features = c("percent.mt"), ncol = 1)
vln_plot3
ggsave("after-qc-percent_mt.png", vln_plot3, dpi = 600)


#Normalization of the data
##The NormalizeData step is basically just ensuring expression values across cells are on a comparable scale. 
##By default, it will divide counts for each gene by the total counts in the cell, multiply that value for each gene by the scale.factor (10,000 by default), and then natural log-transform them.
pc = NormalizeData(object = pc, normalization.method = "LogNormalize", scale.factor = 1e4)
##Finding variable features (2000 highly variable genes that have very high expression in some cells and very low in some)
pc = FindVariableFeatures(pc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

ggsave("highly-variable-features.png", plot2, dpi = 600, bg = "white")

#Scaling the data
pc = ScaleData(pc, verbose = F)

#High dimension reduction (PCA, Heatmap, JackStarwPlot)
pc = RunPCA(pc, npcs = 50, features=VariableFeatures(object = pc), verbose = F)
print(pc[["pca"]], dims = 1:5, nfeatures = 10)

pca_plot <- DimPlot(pc, reduction = "pca", dims = (c(1,2)))
pca_plot 
ggsave("pca_plot.png", pca_plot, dpi = 600, bg = "white")

##plot not saved
DimHeatmap(pc, dims = 1:2, cells = 500, balanced = TRUE)
ggsave("basic_heatmap.png", last_plot(), dpi = 600, bg = "white")


pc = JackStraw(pc, num.replicate = 100)
pc = ScoreJackStraw(pc, dims = 1:20)

options(repr.plot.height = 5, repr.plot.width = 7)
jack_straw <- JackStrawPlot(pc, dims = 1:20)
ggsave("jackstraw_plot.png", last_plot(), dpi = 600, bg = "white")

#DON'T RUN THIS CODE
source('https://raw.githubusercontent.com/constantAmateur/MCVR/master/code.R')
# define a function that does standard processing that is required by the molecularCrossValidation function.
stdProcess = function(dat,...){
  ss = CreateSeuratObject(dat)
  ss = NormalizeData(ss,verbose = F)
  return(ss@assays$RNA@scale.data)
}
#run the function
mcv = molecularCrossValidation(roundToInt(pc@assays$RNA@counts),VariableFeatures(pc),normalisation=stdProcess,features=VariableFeatures(pc),tFracs=c(1,0.9,0.8,0.5),nSplits=5)
#get the optimal number of PCs
max(mcv$mins1sd[[1]])

#Cell clustering based on gene expression data based on nearest-neighbor method implemented in Seurat
pc = FindNeighbors(pc, dims = 1:20, verbose = F)
pc = FindClusters(pc, resolution = 0.6, verbose = F)
pc = RunUMAP(pc, dims = 1:20, verbose = F)

DimPlot(pc, label=T)
ggsave("umap_unannotated.png", last_plot(), dpi = 600, bg = "white")

DimPlot(pc, group.by = "orig.ident", label = F)
ggsave("umap_unannotated_normal_tumor.png", last_plot(), dpi = 600, bg = "white")

pc = RunTSNE(pc, dims = 1:20, verbose = F)
DimPlot(object = pc, reduction = "tsne", label = T)
ggsave("tsne_unannotated.png", last_plot(), dpi = 600, bg = "white")

DimPlot(object = pc, reduction = "tsne", group.by = "orig.ident")
ggsave("tsne_unannotated_normal_tumor.png", last_plot(), dpi = 600, bg = "white")


#optional 
VizDimLoadings(pc, dims = 1:2, reduction = "pca")

#Cell annotation using a supervised method that uses experimentally validated marker database, scType
library(HGNChelper)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = pc[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pc[["RNA"]]@scale.data (default), pc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(pc@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(pc@meta.data[pc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
as.data.frame(cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores))[,2]

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
scores <- sctype_scores[,1:3]
write.csv(scores, "cell_annotations.csv")

pc@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  pc@meta.data$customclassif[pc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(pc_subset, reduction = "umap", label = TRUE, group.by = 'customclassif', repel = T)
ggsave("umap_annotated_normal_tumor_cells_removed.png", last_plot(), dpi = 600, bg = "white")
DimPlot(pc_subset, reduction = "tsne", label = TRUE, group.by = 'customclassif', repel = T)        
ggsave("tsne_annotated_normal_tumor_cells_removed.png", last_plot(), dpi = 600, bg = "white")


#Differential gene expression analysis and subsequent visualizations

levels(pc)
pc <- RenameIdents(object = pc, 
                                  "0" = "CD8+ NKT-like cells",
                                  "1" = "Myeloid Dendritic cells",
                                  "2" = "CD8+ NKT-like cells",
                                  "3" = "Naive CD4+ T cells",
                                  "4" = "Naive CD4+ T cells",
                                  "5" = "Memory CD4+ T cells",
                                  "6" = "Classical Monocytes",
                                  "7" = "Naive B cells",
                      "8" = "CD4+ NKT-like cells",
                      "9" = "CD8+ NKT-like cells",
                      "10" = "Progenitor cells",
                      "11" = "Basophils",
                      "12" = "CD4+ NKT-like cells",
                      "13" = "Plasmacytoid Dendritic cells",
                   "14" = "Cancer cells")

Idents(pc)

# Find differentially expressed features between 
cancer.markers <- FindMarkers(pc, ident.1 = "Cancer cells", ident.2 = "Progenitor cells")
# view results
head(cancer.markers)

cancer.markers <- FindMarkers(, ident.1 = "Cancer cells", ident.2 = "Progenitor cells")
# view results
head(cancer.markers)

#against all other cells
cancer.vs.all <- FindMarkers(pc, ident.1 = "Cancer cells", ident.2 = NULL, only.pos = TRUE)
# view results
head(cancer.vs.all)

##finding all markers, computing altogether 
markers_genes_up <- FindAllMarkers(pc_subset, log2FC.threshold = 0.2, test.use = "wilcox",
    min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE,
    assay = "RNA")

markers_genes_dn <- FindAllMarkers(pc_subset, log2FC.threshold = 0.2, test.use = "wilcox",
                                   min.pct = 0.1, min.diff.pct = 0.2, only.pos = FALSE,
                                   assay = "RNA")

# Filter top 10 genes for each cell type
top10 <- markers_genes_up %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC)

top5 <- markers_genes %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC)

down_top5 <- markers_genes_dn %>%
  group_by(cluster) %>%
  top_n(-5, avg_log2FC)



##gene specific plots
FeaturePlot(object = pc, features = c("GZMB", "GZMH"))
DotPlot(object = pc_subset, features = c("HLA-DRA", "GZMB", "GZMH", "NKG7", "PRF1", "GNLY","HLA-DRB1"), cols = c("lightgrey", "#159895"))
ggsave("genes_expressed_multiple_cells.png", last_plot(), dpi = 600, bg = "white")
DotPlot(object = pc_subset, features = c("TNF", "TNFRSF18", "TNFRSF9", "TNFRSF4"), cols = c("lightgrey", "#159895"))
ggsave("TNF_family.png", last_plot(), dpi = 600, bg = "white")

DotPlot(object = pc_subset, features = c("PLAUR", "CCL5", "CD3D", "CXCL8", "IL32", "TYROBP", "AIF1", "LYZ", "TRAC"), cols = c("lightgrey", "#159895"))
ggsave("common_downregulated_genes.png", last_plot(), dpi = 600, bg = "white")

DotPlot(object = pc_subset, features = c("HLA-DPA1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1"), cols = c("lightgrey", "#159895"))
ggsave("MHC_family.png", last_plot(), dpi = 600, bg = "white")


DotPlot(object = pc, features = c("FOXP3", "IL32"))
FeatureScatter(object = pc, feature1 = "MPO", feature2 = "PC_1")
FeatureScatter(object = pc, feature1 = "MS4A1", feature2 = "CD3D")

VariableFeaturePlot(object = pc)

# Violin and Ridge plots
VlnPlot(object = pc_subset, features = c("HLA-DPA1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5"))
VlnPlot(object = pc, features = c("BRCA1"))
RidgePlot(object = pc, feature = c("HBD"))
RidgePlot(object = pc_subset, feature = c("TNFRSF9"))
ggsave("TNFRSF9.png", last_plot(), dpi = 600, bg = "white")

# Heatmaps
friendly_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99")
DoHeatmap(object = pc_subset, features = top10$gene)
DoHeatmap(subset(pc_subset, downsample = 100), features = top10$gene, size = 3, group.colors = friendly_cols, draw.lines = TRUE, label = FALSE)
ggsave("heatmap_top5_genes.png", last_plot(), dpi = 900)

write.csv(cancer.markers, "cancer.markers.vs.all.csv")



library(CellChat)
#Let's convert the Seurat-based object into cellChat object
cellchat <- createCellChat(object = pc_subset, group.by = "customclassif")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB) #multiple categories available

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use


# subset the expression data of signaling genes only, ignoring the rest of the genes
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multicore", workers = 24) # to utilize multiple processors/threads to speed up the work

#overexpressed genes are being identified
cellchat <- identifyOverExpressedGenes(cellchat)
#overexpressed interactions have to be identified, let's do that first
cellchat <- identifyOverExpressedInteractions(cellchat)


# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#let's project the gene expression data into protein-protein interactions data (gene-gene/protein-protein)
cellchat <- projectData(cellchat, PPI.human)

#calculating probabilities of all the protein-protein interactions
cellchat <- computeCommunProb(cellchat, raw.use = FALSE) #raw.use = false makes sure we use the projected data

#Communication is identified, let's filter bad-quality communications/less number of cells per group
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

#calculating the pathways probability for the cell-cell communications found above
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#Let's do network visualization of the cells and their communications
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
vis1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
png("all_cells_interactions.png", units = "in", res = 600)
print(vis1)
dev.off()
vis2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
png("all_cells_weights_main.png", units = "in", res = 600)
print(vis2)
dev.off()


#Visualization of each cell type against other cell types 1 by 1
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Let's visualize a single pathway
pathways.show <- c("TNF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


# Chord diagram
par(mfrow=c(1,1))
pdf(file =paste0(pathways.show, "-chord.pdf"), width = 20, height =16)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()


# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> 
netAnalysis_contribution(cellchat, signaling = pathways.show)


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.png"), plot=gg, width = 3, height = 2, units = 'in', dpi = 600)
}


#Let's visualize highest number of interactions of VEGF pathway
netAnalysis_contribution(cellchat, signaling = "TNF")

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1:4,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

pdf(file ="cellchat1.pdf", width = 20, height =16)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 8, targets.use = c("Neutrophils"), remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = "CD4+ NKT-like cells", signaling = c("TNF","CXCL", "CSF3"), remove.isolate = FALSE)

pdf(file ="cellchat2.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
dev.off()

pdf(file ="cellchat3.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat, sources.use = c("Neutrophils", "Pro_B Cells", "Platelets"), targets.use = "Cancer cells", legend.pos.x = 15)
dev.off()

plotGeneExpression(cellchat, signaling = "COMPLEMENT")


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = "TNF", width = 8, height = 2.5, font.size = 10)


ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("TNF", "CCL"))
ht

####important
netVisual_bubble(cellchat, sources.use = 1, signaling = c("CCL","CXCL", "CLEC", "MHC-I", "MHC-II"))
ggsave("cd4+_signaling.png", last_plot(), dpi = 900)

netVisual_bubble(cellchat, sources.use = 2, signaling = c("CD99", "CDH1", "CCL","CXCL", "CLEC", "MHC-I", "MHC-II"))
ggsave("cd8+_signaling.png", last_plot(), dpi = 900)

netVisual_bubble(cellchat, sources.use = 1, signaling = c("CD99", "CDH1", "CCL","CXCL", "CLEC", "MHC-I", "MHC-II", "CSF3", "CSF"))
ggsave("naive_cd4_cells_signaling.png", last_plot(), dpi = 900)


netVisual_chord_gene(cellchat, sources.use = 2, signaling = c("CD99", "CDH1", "CCL","CXCL", "CLEC", "MHC-I", "MHC-II", "CSF3", "CSF"),legend.pos.x = 8)


ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = c("CD99", "CDH1", "CCL","CXCL", "CLEC", "MHC-I", "MHC-II", "CSF3", "CSF"))
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = c("CD99", "CDH1", "CCL","CXCL", "CLEC", "MHC-I", "MHC-II", "CSF3", "CSF"))
ht1 + ht2


# List of pathways
pathways <- c("MHC-I", "MHC-II", "CCL", "CLEC", "CXCL", "CD99", "GALECTIN", "ADGRE5", "IL2", "MIF", "ICAM", "IL1", "LCK", "ITGB2", "ANNEXIN", "VISFATIN", "LT", "RESISTIN", "SPP1", "IL16", "CD96", "BAFF", "TNF", "COMPLEMENT", "FN1", "TIGIT", "PVR", "APP", "IL4", "OX40", "CD86", "CSF", "BTLA", "GAS", "NECTIN", "LIGHT", "THBS", "CD80", "CD40", "CD22", "CD45", "VEGF", "GRN", "IL6", "OSM", "VEGI", "FASLG", "CD48", "GP1BA", "CD23", "CSF3", "L1CAM", "PARs", "CDH1", "LIFR", "ALCAM", "CD6", "NPR2", "THY1", "IGF", "PDGF", "SELL", "MPZ", "NOTCH", "CNTN", "HGF", "EPHA", "BMP", "EPHB", "CDH", "ACTIVIN", "FGF", "FLT3", "OCLN", "WNT", "SEMA4")

# Create a directory to save the plots
dir.create("pathway_plots")

# Loop through each pathway and save the gene expression plot
for (pathway in pathways) {
  # Generate the plot
  plotGeneExpression(cellchat, signaling = pathway)
  
  # Set the filename for saving the plot
  filename <- paste0("pathway_plots/", pathway, ".png")
  
  # Save the plot with dpi 600 and white background
  ggsave(filename, dpi = 600, bg = "white")
}


# Assuming your Seurat object is named 'seurat_obj'
# Remove cells belonging to specific groups based on a marker gene

groups_to_remove <- c("Basophils", "Plasmacytoid Dendritic cells", "Cancer cells", "Progenitor cells")  # Replace with your actual group names

subset_pc <- subset(pc, subset = !(Idents(seurat_obj) %in% groups_to_remove))

pc_subset <- pc[, !(pc@meta.data$customclassif %in% groups_to_remove)]

