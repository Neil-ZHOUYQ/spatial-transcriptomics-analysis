# This script load data to create seurat objects, investigate the spots , do quality control and filters spots.





# Libraries
snapshot <- "/path_to_R_snapshot/"
.libPaths(snapshot)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

# Set working directory
setwd("/your_working_directory/")

########## 1. Load the raw dataset ##########
samples_healthy <- c(25, 62)
samples_ipf <- c(82, 7)
samples <- c(samples_healthy, samples_ipf)

# Create spatial trans objects
Spatialtrans_list <- list()


for (sample in samples) {
    cat(paste0("Loading spatial dataset for sample ", sample, "...\n"))
    data_dir <- paste0("/path_to_STdata/", sample, "/outs")  
    Spatialtrans <- Load10X_Spatial(
        data.dir = data_dir,
        filename = "filtered_feature_bc_matrix.h5",
        assay = "Spatial",
        slice = paste0("slice", sample),
        filter.matrix = TRUE,
        to.upper = FALSE)
    Spatialtrans_list[[paste0("Lung_", sample)]] <- Spatialtrans
}

# p220 samples
samples_p220 <- c(23, 57)
Spatialtrans_p220 <- list()
for (sample in samples_p220) {
    cat(paste0("Loading spatial dataset for sample ", sample, "...\n"))
    data_dir <- paste0("/path_to_spaceranger_output/", sample, "/outs")  
    Spatialtrans <- Load10X_Spatial(
        data.dir = data_dir,
        filename = "filtered_feature_bc_matrix.h5",
        assay = "Spatial",
        slice = paste0("slice", sample),
        filter.matrix = TRUE,
        to.upper = FALSE)
    Spatialtrans_p220[[paste0("Lung_", sample)]] <- Spatialtrans
}

# The visium data from 10x consists of the following data types:
# 1. A spot by gene expression matrix
# 2. An image of the tissue slice (obtained from H&E staining during data acquisition)
# 3. Scaling factors that relate the original high resolution image to the lower resolution image used here for visualization.

# All samples
samples <- c(samples, samples_p220)
Spatialtrans_list <- c(Spatialtrans_list, Spatialtrans_p220)      # seurat objects
Spatialtrans_list

########## 2. Data Preprocessing ##########
# The initial preprocessing steps that we perform on the spot by gene expression data are similar to a typical scRNA-seq experiment. We first need to normalize the data in order to account for variance in sequencing depth across data points. We note that the variance in molecular counts / spot can be substantial for spatial datasets, particularly if there are differences in cell density across the tissue. We see substantial heterogeneity here, which requires effective normalization.
Spatialtrans_list_filtered <- list()
NumofSpots <- matrix(0, length(samples), 2)
rownames(NumofSpots) <- paste0("Lung_", samples)
colnames(NumofSpots) <- c("Prefiltering", "Postfiltering")
for (sample in samples) {
  cat("QC for sample", sample, "\n")
  spatialtrans_obj <- Spatialtrans_list[[paste0("Lung_", sample)]]

  # Wrap plots pre filtering
  p1 <- VlnPlot(spatialtrans_obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  p2 <- SpatialFeaturePlot(spatialtrans_obj, features = "nCount_Spatial") + theme(legend.position = "right")      
  p <- wrap_plots(p1, p2)
  png(paste0("Figures/wrap_plots_", sample, "_raw.png"), res = 300, height = 1500, width = 3000)
  print(p)  
  dev.off()
  
  # QC: Mitochondria percentage
  spatialtrans_obj <- PercentageFeatureSet(spatialtrans_obj, "^MT-", col.name = "percent_mito")
  NumofSpots[paste0("Lung_", sample), "Prefiltering"] <- ncol(spatialtrans_obj)
  cat("Number of spots for sample", sample, "prefiltering:", ncol(spatialtrans_obj), "\n")

  # Histograms of QC metrics
  png(paste0("Figures/Histograms_QCmetrics_", sample, ".png"), res = 300, height = 1200, width = 3600)
  par(mfrow = c(1, 3))
  hist(log10(spatialtrans_obj$nCount_Spatial), 20, xlab = "log10_UMIs per spot", main = "nCount_Spatial", col = "lightblue", border = F)
  hist(log10(spatialtrans_obj$nFeature_Spatial), 20, xlab = "log10_Genes per spot", main = "nFeature_Spatial", col = "lightblue", border = F)
  hist(spatialtrans_obj$percent_mito, 20, xlab = "Percent mito UMIs", main = "percent_mito", col = "lightblue", border = F)
  dev.off()

  p <- ggplot(spatialtrans_obj@meta.data, aes(nCount_Spatial, nFeature_Spatial)) + geom_point(size=2) + scale_x_log10() + scale_y_log10()
  png(paste0("Figures/ScatterPlot_nCountsVSnFeatures_", sample, ".png"), res = 300, height = 1500, width = 1500)
  print(p)
  dev.off()
  
  # QC violin plot pre-filtering
  p <- VlnPlot(
    object = spatialtrans_obj[, spatialtrans_obj$nCount_Spatial > 0 & spatialtrans_obj$nFeature_Spatial > 0],
    features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"),
    ncol = 3,
    pt.size = 0.01
  )
  png(paste0("Figures/Vlnplot_", sample, "_prefiltering.png"), res = 300, height = 1500, width = 3000)
  print(p)
  dev.off()

  # QC violin plot after filtering: Use nFeatures larger or equal to 50 and mitochodria percentage smaller than 20% as the threshold
  spatialtrans_obj <- spatialtrans_obj[, spatialtrans_obj$nFeature_Spatial >= 50 & spatialtrans_obj$percent_mito < 20]                      # filtering step: keep spots satisfies requirements
  
  
  
  # QC violin plot after filtering: Use nFeatures larger or equal to 50 and mitochodria percentage smaller than 20% as the threshold
  p <- VlnPlot(
    object = spatialtrans_obj,
    features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"),
    ncol = 3,
    pt.size = 0.01
  )
  png(paste0("Figures/Vlnplot_", sample, "_postfiltering.png"), res = 300, height = 1500, width = 3000)
  print(p)
  dev.off()

  # Wrap plots post filtering
  p1 <- VlnPlot(spatialtrans_obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  p2 <- SpatialFeaturePlot(spatialtrans_obj, features = "nCount_Spatial") + theme(legend.position = "right")
  p <- wrap_plots(p1, p2)
  png(paste0("Figures/wrap_plots_", sample, "_filtered.png"), res = 300, height = 1500, width = 3000)
  print(p)  
  dev.off()

  # Update Spatialtrans_list
  NumofSpots[paste0("Lung_", sample), "Postfiltering"] <- ncol(spatialtrans_obj)
  cat("Number of spots for sample", sample, "postfiltering:", ncol(spatialtrans_obj), "\n")
  Spatialtrans_list_filtered[[paste0("Lung_", sample)]] <- spatialtrans_obj
}
NumofSpots
Spatialtrans_list_filtered
saveRDS(Spatialtrans_list_filtered, file = "IntermediateResults/Spatialtrans_list_filtered.rds")

