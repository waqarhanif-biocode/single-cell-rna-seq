install.packages("devtools")
devtools::install_github("sqjin/CellChat")
install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
#python-package: pip install umap-learn

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#Read the scRNA-seq dataset:
#to read 10x data specifically
tumor_data = "filtered" #directory where you have the gene-count matrix
tumor_data = Read10X(data.dir = tumor_data)  # Seurat function to read in 10x count data
dim(tumor_data)

tumor = CreateSeuratObject(tumor_data, project = "CML") #give your project a name
tumor@meta.data$disease = "CML_tumor" #add a metadata related to the phenotype


#Let's convert the Seurat-based object into cellChat object
cellchat <- createCellChat(object = tumor, group.by = "customclassif")

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
future::plan("multiprocess", workers = 4) # to utilize multiple processors/threads to speed up the work

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
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Visualization of each cell type against other cell types 1 by 1
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Let's visualize a single pathway
pathways.show <- c("CSF3") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


# Chord diagram
par(mfrow=c(1,1))
pdf(file ="cellchat.pdf", width = 20, height =16)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()


# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
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
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


#Let's visualize highest number of interactions of VEGF pathway
netAnalysis_contribution(cellchat, signaling = "MIF")

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
netVisual_bubble(cellchat, sources.use = "Cancer cells", targets.use = c("Neutrophils", "Pro_B Cells", "Platelets", "Unknown"), signaling = c("VEGF","CXCL", "CD34", "CSF3"), remove.isolate = FALSE)

pdf(file ="cellchat2.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
dev.off()

pdf(file ="cellchat3.pdf", width = 20, height =16)
netVisual_chord_gene(cellchat, sources.use = c("Neutrophils", "Pro_B Cells", "Platelets"), targets.use = "Cancer cells", legend.pos.x = 15)
dev.off()

plotGeneExpression(cellchat, signaling = "MIF")


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = "MIF", width = 8, height = 2.5, font.size = 10)


ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CSF3", "MIF"))
ht
