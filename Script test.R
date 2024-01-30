#------------------------------------------------------------------------------#
# Pipeline Single cell analysis Panc8
#------------------------------------------------------------------------------#
# from: https://www.singlecellcourse.org/, 
# https://satijalab.org/seurat/articles/integration_mapping.html


# 1. Install and Load required packages.========================================

# # Devtools.---------------------------------------------------------------------
# install.packages('devtools')
# library(devtools)
# 
# # Seuratdata.-------------------------------------------------------------------
# devtools::install_github('satijalab/seurat-data')
# 
# # Seurat wrappers.--------------------------------------------------------------
# devtools::install_github('satijalab/seurat-wrappers')

# Load packages.----------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(ggplot2)

# library(SeuratWrappers)
# library(patchwork)

# 2. Load Data.=================================================================
# Load Panc8 data.--------------------------------------------------------------

# InstallData('panc8')

panc8 <- LoadData("panc8")
panc8 <- UpdateSeuratObject(object = panc8)   # Update the object.

head(panc8[[]])           # Seurat metadata.
table(panc8$tech)         # All tech used.

# 3. Select data to analyze.====================================================

pancreas.ref <- subset(panc8, tech %in% c("celseq2", "smartseq2"))
class(pancreas.ref[["RNA"]])                  # Change from Assay to Assay5.
pancreas.ref[["RNA"]] <- as(object = pancreas.ref[["RNA"]], Class = "Assay5")

pancreas.ref[["RNA"]] <- split(pancreas.ref[["RNA"]], f = pancreas.ref$tech)

# 4. Pre-process dataset (without integration).=================================

pancreas.ref <- NormalizeData(pancreas.ref)
pancreas.ref <- FindVariableFeatures(pancreas.ref)
pancreas.ref <- ScaleData(pancreas.ref)
pancreas.ref <- RunPCA(pancreas.ref)
  DimPlot(pancreas.ref)  # Plot PCA.
pancreas.ref <- FindNeighbors(pancreas.ref, dims = 1:30)
pancreas.ref <- FindClusters(pancreas.ref)

# 5. Visualize UMAP before integration.=========================================
pancreas.ref <- RunUMAP(pancreas.ref, dims = 1:30)
DimPlot(pancreas.ref, group.by = "tech")

# 6. Integrate Layers.==========================================================
pancreas.ref <- IntegrateLayers(object = pancreas.ref, method = CCAIntegration, orig.reduction = "pca",
                                new.reduction = "integrated.cca", verbose = FALSE)
pancreas.ref <- FindNeighbors(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
pancreas.ref <- FindClusters(pancreas.ref)

# 7. Visualize UMAP after integration.==========================================
pancreas.ref <- RunUMAP(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
DimPlot(pancreas.ref, group.by = c("tech", "celltype"))

# 8. Cell type classification using an integrated reference.====================
# select two technologies for the query datasets
pancreas.query <- subset(panc8, tech %in% c("fluidigmc1", "celseq"))
pancreas.query <- NormalizeData(pancreas.query)
pancreas.anchors <- FindTransferAnchors(reference = pancreas.ref, 
                                        query = pancreas.query, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, 
                            refdata = pancreas.ref$celltype, dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

# Compare to the original label annotations of full integrated analysis.--------
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)        # 96% of cells labeled correctly.

# Verify labels.----------------------------------------------------------------
table(pancreas.query$predicted.id)

# 9. Check expression of genes in all populations.==============================

# Violin Plot.------------------------------------------------------------------

VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), 
        group.by = "predicted.id")

# Genes specifically expressed in islet cells of pancreas
VlnPlot(pancreas.query,features = c("INS","GCG","IAPP"), 
        group.by = "predicted.id")
# Genes specifically expressed in exocrine glandular cells of pancreas
VlnPlot(pancreas.query, features = c("AMY2A","CELA3A","CPA1"), 
        group.by = "predicted.id")
# Genes specifically expressed in ductal cells of pancreas
VlnPlot(pancreas.query, features = c("SCTR","CFTR","DCDC2"), 
        group.by = "predicted.id")


# End of script.
################################################################################