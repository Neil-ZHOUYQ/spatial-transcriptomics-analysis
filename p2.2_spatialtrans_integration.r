# This script do normalization, romove batch effect, PCA, select dims and resolution, and view dimplot of spots.




# Libraries
snapshot <- "/path_to_R_snapshot/"
.libPaths(snapshot)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

# Set working directory
setwd("/your_working_directory/")

# Load the filtered spatial trans object
Spatialtrans_list_filtered <- readRDS("IntermediateResults/Spatialtrans_list_filtered.rds")
samples <- names(Spatialtrans_list_filtered)

##### 1. SCTransform normalization
cat("SCTranform normalization for all samples...\n")
for (sample in samples) {
    cat("SCTranform for", sample, "...\n")
    Spatialtrans_list_filtered[[sample]] <- SCTransform(Spatialtrans_list_filtered[[sample]], assay = "Spatial", verbose = FALSE)
}

##### 2. Sample integration
# Rename the cell names to enforce unique cell names
for (sample in samples) {
  Spatialtrans_list_filtered[[sample]]$orig.ident <- sample
  Spatialtrans_list_filtered[[sample]] <- RenameCells(Spatialtrans_list_filtered[[sample]], add.cell.id = sample)   # make sure each spot has unique ID
}

# Select features
cat("Select features...\n")
features <- SelectIntegrationFeatures(object.list = Spatialtrans_list_filtered, nfeatures = 3000, verbose = FALSE)  # select 3000 HVGs

# Prep for integration
cat("Prepare for integration...\n")
options(future.globals.maxSize = 2 * 1024^3)
PrepSCT <- PrepSCTIntegration(object.list = Spatialtrans_list_filtered, anchor.features = features, verbose = FALSE)

# Find anchors, which are spots representing the same biological items, core step of integration
cat("Find anchors...\n")
integ.anchors <- FindIntegrationAnchors(object.list = PrepSCT, normalization.method = "SCT", verbose = FALSE, anchor.features = features) 

# Integrate
cat("SCT integration...\n")
Spatialtrans.integrated <- IntegrateData(anchorset = integ.anchors, normalization.method = "SCT", verbose = FALSE)   # one seurat object containing data from 6 samples, with batch effect corrected

# Run PCA and generate elbow plots
pca_npcs <- 50
cat("Run PCA and generate elbow plots...\n")
Spatialtrans.integrated <- RunPCA(Spatialtrans.integrated, npcs = pca_npcs, verbose = FALSE)
png("Figures/Elbowplot_integrated.png", res = 300, height = 1500, width = 1500)
ElbowPlot(Spatialtrans.integrated, ndims = pca_npcs)         # to decide the PCs that we use
dev.off()
saveRDS(Spatialtrans.integrated, file = "IntermediateResults/Spatialtrans.integrated.rds")

##### 3. Determine the dimension and resolution for cell annotation
# Run permutations of dim and res to determine the optimal dimension and resolution
dimensions = c(15, 20, 25)  
resolutions = c(0.5, 0.8, 1, 1.2, 1.5, 2)
dir.create("Results/rds/")
for(dim_num in dimensions){
    Spatialtrans.integrated.perm <- Spatialtrans.integrated
    for(res_num in resolutions){
        cat("Processing dimension", dim_num, "resolution", res_num, "...\n")
        Spatialtrans.integrated.perm <- FindNeighbors(Spatialtrans.integrated.perm, reduction = "pca", dims = 1:dim_num)
	      Spatialtrans.integrated.perm <- FindClusters(Spatialtrans.integrated.perm, resolution = res_num, random.seed = 0)               # clustering
	      Spatialtrans.integrated.perm <- RunUMAP(Spatialtrans.integrated.perm,  reduction = "pca", dims = 1:dim_num)
	      Spatialtrans.integrated.perm <- RunTSNE(Spatialtrans.integrated.perm, reduction = "pca", dims = 1:dim_num, seed.use = 1)
    }

    # change to factor
    rel.cols <- grep("integrated_snn",colnames(Spatialtrans.integrated.perm@meta.data))
    for(my.col in rel.cols){
        Spatialtrans.integrated.perm@meta.data[,my.col] <- as.character(Spatialtrans.integrated.perm@meta.data[,my.col])
    }

    # Saving RDS file
    rds_file <- paste0("Results/rds/spatialtrans_d",dim_num,".rds")
    saveRDS(Spatialtrans.integrated.perm, file = rds_file)

    # Saving H5AD
    SaveH5Seurat(Spatialtrans.integrated.perm, overwrite = T, filename = paste0("Results/rds/spatialtrans_d",dim_num,".h5Seurat"))
    Convert(paste0("Results/rds/spatialtrans_d",dim_num,".h5Seurat"), overwrite = T, dest = "h5ad")             # convert to .h5ad, which can be analyzed by python scanpy
}

# Seleted dimension and resolution
dim_num = 20
res_num = 0.5
print(paste0("Selected optimal dimension for downstream processing d=",dim_num,"; r=",res_num))

# Read the seurat obj
Spatialtrans.integrated <- readRDS(paste0("Results/rds/spatialtrans_d", dim_num, ".rds"))
Idents(Spatialtrans.integrated) <- Spatialtrans.integrated@meta.data[, paste0("integrated_snn_res.", res_num)] # use cluster results of selected resolution as the active ident 
Spatialtrans.integrated <- AddMetaData(Spatialtrans.integrated, metadata = Spatialtrans.integrated$orig.ident, col.name = "group")
Spatialtrans.integrated$group <- plyr::mapvalues(Spatialtrans.integrated$group, from = c("Lung_23", "Lung_25", "Lung_62", "Lung_7", "Lung_57", "Lung_82"), to = c("Healthy", "Healthy", "Healthy", "IPF", "IPF", "IPF")) # transform group names from sample names to "Healthy"&"IPF"

# Sort the levels
levels(Spatialtrans.integrated) <- seq(0, length(levels(Spatialtrans.integrated))-1)

# UMAP plot
print("Generating Initial UMAP plot")
p <- DimPlot(Spatialtrans.integrated, reduction = "umap", label = TRUE, repel = TRUE)
png(paste0("Figures/UMAP_spatialtrans_d", dim_num, "_r", res_num, ".png"), res = 300, width = 1750, height = 1750)
print(p)
dev.off()
p <- DimPlot(Spatialtrans.integrated, reduction = "umap", label = TRUE, repel=TRUE, split.by = "group")
png(paste0("Figures/UMAP_spatialtrans_d", dim_num, "_r", res_num, "_splitbyconds.png"), res = 300, width = 3000, height = 1500)
print(p)
dev.off()
p <- DimPlot(Spatialtrans.integrated, reduction = "umap", label = TRUE, repel=TRUE, split.by = "orig.ident", ncol = 3)
png(paste0("Figures/UMAP_spatialtrans_d", dim_num, "_r", res_num, "_splitbysamples.png"), res = 300, width = 3500, height = 2250)
print(p)
dev.off()

# Spatial plot, core of spacial trancriptomic
# Healthy
Spatialtrans.integrated_subset <- Spatialtrans.integrated[, Spatialtrans.integrated$group == "Healthy"]      
Spatialtrans.integrated_subset@images[["slice7"]] = NULL
Spatialtrans.integrated_subset@images[["slice57"]] = NULL
Spatialtrans.integrated_subset@images[["slice82"]] = NULL
p <- SpatialDimPlot(Spatialtrans.integrated_subset, ncol = 3)
png(paste0("Figures/spatialdimplot_d", dim_num, "_r", res_num, "_healthy.png"), res = 300, width = 4500, height = 2000)
print(p)
dev.off()
# IPF
Spatialtrans.integrated_subset <- Spatialtrans.integrated[, Spatialtrans.integrated$group == "IPF"]
Spatialtrans.integrated_subset@images[["slice23"]] = NULL
Spatialtrans.integrated_subset@images[["slice25"]] = NULL
Spatialtrans.integrated_subset@images[["slice62"]] = NULL
p <- SpatialDimPlot(Spatialtrans.integrated_subset, ncol = 3)
png(paste0("Figures/spatialdimplot_d", dim_num, "_r", res_num, "_ipf.png"), res = 300, width = 4500, height = 2000)
print(p)
dev.off()

# Create spatial embedding
all_coords <- list()
for (slice in names(Spatialtrans.integrated@images)) {
  coords <- Spatialtrans.integrated@images[[slice]]@coordinates[, c("row", "col")]
  all_coords[[slice]] <- as.matrix(coords)
}
spatial_coords <- do.call(rbind, all_coords)
spatial_coords <- spatial_coords[colnames(Spatialtrans.integrated), ] #order
colnames(spatial_coords) <- paste0("SPATIAL_", 1:2)
Spatialtrans.integrated[["spatial"]] <- CreateDimReducObject(embeddings = spatial_coords, key = "SPATIAL_", assay = DefaultAssay(Spatialtrans.integrated))  # can use DimPlot(Spatialtrans.integrated, reduction = "spatial")

# Save the current integrated Seurat version
saveRDS(Spatialtrans.integrated, paste0("Results/Spatial_Integrated_d", dim_num, "_r", res_num, ".rds"))

