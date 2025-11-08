# use SPOtlight to do deconvolution, annotation


# Libraries
snapshot <- "/path_to_R_snapshot/"
.libPaths(snapshot)
library(Seurat)
library(SeuratDisk)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(ComplexHeatmap)
library(scProportionTest)

# Set working directory
setwd("/your_working_directory/")

# Load the integrated spatial transcriptomic dataset
cat("Loading the integrated spatial dataset and lung reference.\n")
dim_num = 20
res_num = 0.5
Spatial_Integrated <- readRDS(paste0("Results/Spatial_Integrated_d", dim_num, "_r", res_num, ".rds"))

cat("Loading the scRNAseq reference dataset...\n")
# Reference data: LungCellAtlas (processed for p290)
refdata <- "LungCellAtlas"
lung_reference_full  <- readRDS("/path_to_scref/scref.rds")

# Use 'ann_level_2' or 'ann_level_3' for cell type annotation, how precise the annotation is 
celltype2use = "ann_level_3"
print(paste("Pre-Filtering Table for", celltype2use))
print(table(lung_reference_full@meta.data[, celltype2use]))

# Remove NAs, None, Unkonwns and celltypes less than 100 cells
idx2remove <- grep(TRUE,is.na(lung_reference_full$ann_level_3))
lung_reference_full <- lung_reference_full[, -idx2remove]
idx2remove <- which(lung_reference_full@meta.data[, celltype2use] %in% c("None","Unknown","Lymphatic EC proliferating"))
lung_reference_full <- lung_reference_full[, -idx2remove]
print(paste("Post-Filtering Table for", celltype2use))
print(table(lung_reference_full@meta.data[, celltype2use]))

# Create SingleCellExperiment object for reference single cell data and SpatialExperiment object for spatial transcriptomic data
cat("Create SingleCellExperiment object for reference single cell data and SpatialExperiment object for spatial transcriptomic data.\n")
lung_reference_full.sce <- as.SingleCellExperiment(lung_reference_full)
DefaultAssay(Spatial_Integrated) <- "integrated"
Spatial_Integrated.spe <- SpatialExperiment(
  assay = Spatial_Integrated@assays$Spatial@counts, 
  colData = Spatial_Integrated@reductions$spatial@cell.embeddings, 
  spatialCoordsNames = colnames(Spatial_Integrated@reductions$spatial@cell.embeddings)
)
assayNames(Spatial_Integrated.spe) <- "counts"

# Normalizing counts
lung_reference_full.sce <- logNormCounts(lung_reference_full.sce)

# Variance modelling
# We aim to identify highly variable genes that drive biological heterogeneity. By feeding these genes to the model we improve the resolution of the biological structure and reduce the technical noise.
# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(lung_reference_full.sce))
dec <- modelGeneVar(lung_reference_full.sce, subset.row = genes)    # find high variance genes
png(paste("Figures/SPOTlight_Variancemodelling_", refdata, "_", celltype2use, ".png", sep = ""), res = 300, height = 1600, width = 1600)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
dev.off()

# Get the top 3000 genes
hvg <- getTopHVGs(dec, n = 3000)

# Extract the reference data labeling
colLabels(lung_reference_full.sce) <- colData(lung_reference_full.sce)[[celltype2use]]

# Compute marker genes of each cell type
mgs <- scoreMarkers(lung_reference_full.sce, subset.row = genes)

# Then we want to keep only those genes that are relevant for each cell identity:
# Examples include avgLogFC, AUC, pct.expressed, p-value…
mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.6, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)   # a filtered marker genes list

# Deconvolution
res <- SPOTlight(
    x = lung_reference_full.sce,
    y = Spatial_Integrated.spe,
    groups = as.character(lung_reference_full.sce[[celltype2use]]),
    mgs = mgs_df,
    hvg = hvg,                                      # V = W * H
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")

# Save the results after each run
saveRDS(res, file = paste("Results/SPOTlight_res_", refdata, "_", celltype2use, ".rds", sep = ""))

# Extract deconvolution matrixs
res <- readRDS(paste("Results/SPOTlight_res_", refdata, "_", celltype2use, ".rds", sep = ""))
mat <- res$mat

# Assign cell type labels to each spot
pred_clabels <- data.frame(celltype = colnames(mat)[apply(mat, 1, which.max)], pred = apply(mat, 1, max))  # each row is a spot
pred_clabels <- cbind(cellid = rownames(pred_clabels), pred_clabels)                                                                        # spot barcode
png(paste("Figures/Histogram_SPOTlight_pred_", refdata, "_", celltype2use, ".png", sep = ""), res = 300, height = 1600, width = 1600)
hist(pred_clabels$pred, 20, las = 1, main = "Highest Prediction per Spot", xlab = "prediction")
dev.off()

# Cell type annotation
Spatial_Integrated[[paste(refdata, "_", celltype2use, "_SPOTlight", sep = "")]] <- pred_clabels$celltype
cell_type_frequencies <- sort(table(Spatial_Integrated[[paste(refdata, "_", celltype2use, "_SPOTlight", sep = "")]])/ncol(Spatial_Integrated), decreasing = T)
Idents(Spatial_Integrated) <- Spatial_Integrated[[paste(refdata, "_", celltype2use, "_SPOTlight", sep = "")]]
saveRDS(Spatial_Integrated, file = paste0("Results/Spatial_Integrated_", refdata, "_", celltype2use, "_Annotated.rds"))

# Visualize cell types with at least 1% spots
# over_1_percent_types <- names(cell_type_frequencies[cell_type_frequencies > 0.01])
# Spatial_Integrated_subset <- Spatial_Integrated[, as.matrix(Spatial_Integrated[[paste(refdata, "_", celltype2use, "_SPOTlight", sep = "")]]) %in% over_1_percent_types]
Spatial_Integrated <- readRDS(paste0("Results/Spatial_Integrated_", refdata, "_", celltype2use, "_Annotated.rds"))
Spatial_Integrated_subset <- Spatial_Integrated

# UMAP with cell type annotation
# p <- DimPlot(Spatial_Integrated, reduction = "umap", group.by = "ident", cells = WhichCells(Spatial_Integrated, idents = over_1_percent_types))
p <- DimPlot(Spatial_Integrated_subset, reduction = "umap", group.by = "ident")
png(paste("Figures/UMAP_integrated_SPOTlight_", refdata, "_", celltype2use, ".png", sep = ""), res = 300, height = 1500, width = 3000)
print(p)
dev.off()

# Spatial plot with cell type annotation
# Healthy
Spatialtrans.integrated_subset <- Spatial_Integrated_subset[, Spatial_Integrated_subset$group == "Healthy"]
Spatialtrans.integrated_subset@images[["slice7"]] = NULL
Spatialtrans.integrated_subset@images[["slice57"]] = NULL
Spatialtrans.integrated_subset@images[["slice82"]] = NULL
p <- SpatialDimPlot(Spatialtrans.integrated_subset, ncol = 3)
png(paste0("Figures/spatialdimplot_d", dim_num, "_r", res_num, "_healthy_annotated.png"), res = 300, width = 6500, height = 2000)
print(p)
dev.off()
# IPF
Spatialtrans.integrated_subset <- Spatial_Integrated_subset[, Spatial_Integrated_subset$group == "IPF"]
Spatialtrans.integrated_subset@images[["slice23"]] = NULL
Spatialtrans.integrated_subset@images[["slice25"]] = NULL
Spatialtrans.integrated_subset@images[["slice62"]] = NULL
p <- SpatialDimPlot(Spatialtrans.integrated_subset, ncol = 3)
png(paste0("Figures/spatialdimplot_d", dim_num, "_r", res_num, "_ipf_annotated.png"), res = 300, width = 6500, height = 2000)
print(p)
dev.off()

# Plot the cell type across all clusters distribution
CP <- as.matrix(table(Idents(Spatial_Integrated), Spatial_Integrated@meta.data[, paste0("integrated_snn_res.", res_num)]))
colnames(CP) <- paste0("Cluster_", colnames(CP))
CP <- t(CP)
col_fun <-  circlize::colorRamp2(c(0, mean(CP)), c("gray99", "red"))
hm <- Heatmap(CP, col = col_fun,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = '#',
        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.0f", CP[i, j]), x, y, gp = gpar(fontsize = 10))})
png(paste0("Figures/CCDistr_counts_", refdata, "_", celltype2use, ".png"), width = 3000, height = 2500, res = 300)
draw(hm)
dev.off()
write.csv(CP, file = paste0("Results/CCDistr_counts_", refdata, "_", celltype2use, ".csv"), quote = F)

# Extract the color
color_map <- p$scales$scales[[1]]$palette(1)

# Extract NMF model fit
mod <- res$NMF

# Topic profiles
png(paste("Figures/SPOTlight_TopicProfiles_", refdata, "_", celltype2use, ".png", sep = ""), res = 300, height = 1500, width = 1500)
plotTopicProfiles(
    x = mod,
    y = lung_reference_full.sce[[celltype2use]],
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)
dev.off()

# Topic profiles by cells
png(paste("Figures/SPOTlight_TopicProfiles_bycells_", refdata, "_", celltype2use, ".png", sep = ""), res = 300, height = 3500, width = 2500)
plotTopicProfiles(
    x = mod,
    y = lung_reference_full.sce[[celltype2use]],
    facet = TRUE,
    min_prop = 0.01,
    ncol = 6)
dev.off()

# Genes learned for each topic
library(NMF)
sign <- mod$w
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)
topgenes4topics <- lapply(seq_len(ncol(sign)), function(i) {
    column_data <- sign[, i]
    top_genes <- head(order(column_data, decreasing = TRUE), 5)
    top_genes <- rownames(sign)[top_genes]
    })
names(topgenes4topics) <- colnames(sign)
saveRDS(topgenes4topics, file = paste0("IntermediateResults/topgenes4topics_", refdata, "_", celltype2use, ".rds"))

# Dotplot for the top 5 genes of each cell type
DefaultAssay(Spatial_Integrated) <- "SCT"
cellidsbytype <- split(colnames(Spatial_Integrated), Spatial_Integrated[[paste(refdata, "_", celltype2use, "_SPOTlight", sep = "")]])
result_list <- lapply(cellidsbytype, function(cols) {
  row_means <- rowMeans(Spatial_Integrated[, cols])
  return(row_means)
})
geneexprbycelltype <- do.call(cbind, result_list)
topgenes4topics <- lapply(topgenes4topics, function(vec){vec <- intersect(vec, rownames(geneexprbycelltype))})
genes2remove <- names(which(table(unlist(topgenes4topics))>1))
topgenes4topics <- lapply(topgenes4topics, function(vec){vec <- setdiff(vec, genes2remove)})
topgeneexprbycelltype <- lapply(topgenes4topics, function(rows){col_means <- colMeans(geneexprbycelltype[rows, ])})
topgeneexprbycelltype <- do.call(rbind, topgeneexprbycelltype)

# topgeneexprbycelltype <- topgeneexprbycelltype[, levels(Idents(Spatial_Integrated))]
# max_index <- rownames(topgeneexprbycelltype)[apply(topgeneexprbycelltype, 2, function(col) which.max(col))]
# max_index <- c(unique(max_index), setdiff(rownames(topgeneexprbycelltype), max_index))
max_index <- data.frame(cbind(colnames(topgeneexprbycelltype)[apply(topgeneexprbycelltype, 1, function(row) which.max(row))], apply(topgeneexprbycelltype, 1, function(row) max(row))))
colnames(max_index) <- c("celltype", "meanexpr")
max_index <- max_index[order(max_index$meanexpr, decreasing = TRUE), ]
# max_index <- c(unique(max_index), setdiff(rownames(topgeneexprbycelltype), max_index))
# topgeneexprbycelltype <- topgeneexprbycelltype[max_index, ]
topic_order <- rownames(max_index)
celltype_order <- c(unique(max_index$celltype), setdiff(levels(Idents(Spatial_Integrated)), max_index$celltype))
Spatial_Integrated.dp <- Spatial_Integrated
Idents(Spatial_Integrated.dp) <- factor(Idents(Spatial_Integrated.dp), levels = celltype_order)
png(paste("Figures/DP_SPOTlight_", refdata, "_", celltype2use, ".png", sep = ""), width = 5000, height = 1250, res = 300)
DotPlot(Spatial_Integrated.dp, features = topgenes4topics, dot.min = 0.01, col.min = 0, dot.scale = 4) +
theme_light() + 
theme(panel.grid = element_blank(), 
panel.border = element_blank()) +
xlab('Marker Genes') +
ylab('Cell Types') +
scale_color_gradient2(low = rgb(1,1,1,0), mid = rgb(1,1,1,0)) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))
dev.off()

# Divide the integrated SpatialExperiment object into two samples
# Add the image to the SpatialExperiment object
Spatial_Integrated.spe
data_dir <- c(
    Lung_23 = "/path_to_23/outs",
    Lung_57 = "/path_to_57/outs",
    Lung_25 = "/path_to_25/outs",
    Lung_62 = "/path_to_62/outs",
    Lung_82 = "/path_to_82/outs",
    Lung_7 = "/path_to_7/outs"
)

samples <- unique(Spatial_Integrated@meta.data$orig.ident)
Spatial_Integrated.spe_bysamples <- list()
for (sample in samples) {
    spots <- colnames(Spatial_Integrated)[which(Spatial_Integrated@meta.data$orig.ident == sample)]
    Spatial_Integrated.i <- Spatial_Integrated.spe[, spots]
    Spatial_Integrated.i <- addImg(Spatial_Integrated.i, 
        sample_id = "sample01", 
        image_id = sample,
        imageSource = paste(data_dir[sample], "/spatial/detected_tissue_image.jpg", sep = ""), 
        scaleFactor = NA_real_, 
        load = TRUE)
    
    # Update Spatial_Integrated.spe_bysamples
    Spatial_Integrated.spe_bysamples[[sample]] <- Spatial_Integrated.i
}

# Spatial correlation and Co-localization plots
for (sample in samples) {
    spots <- colnames(Spatial_Integrated)[which(Spatial_Integrated@meta.data$orig.ident == sample)]
    mat.i <- mat[spots, ]

    # Spatial Correlation Matrix
    p <- plotCorrelationMatrix(mat.i) + ggplot2::ggtitle(sample)
    png(paste("Figures/SPOTlight_spatialcorr_", refdata, "_", celltype2use, "_", sample, ".png", sep = ""), res = 300, height = 2000, width = 2000)
    print(p)
    dev.off()

    # Co-localization
    p <- plotInteractions(mat.i, which = "heatmap", metric = "prop") + ggplot2::ggtitle(sample)
    png(paste("Figures/SPOTlight_colocalization_", refdata, "_", celltype2use, "_", sample, ".png", sep = ""), res = 300, height = 2000, width = 2000)
    print(p)
    dev.off()
}

# Install the regular version for spatial pie chart plot
detach("package:SPOTlight")
snapshot <- "/path_to_snapshot/"
.libPaths(snapshot)
# BiocManager::install("SPOTlight", lib=snapshot)
library(SPOTlight)

# Visualize the cell type proportions as sections of a pie chart for each spot
pal <- color_map

# All samples
for (sample in samples) {
    spots <- colnames(Spatial_Integrated)[which(Spatial_Integrated@meta.data$orig.ident == sample)]
    mat.i <- mat[spots, ]
    mat.i <- mat.i[, names(pal)]
    mat.i[mat.i < 0.1] <- 0
    
    Spatial_Integrated.i <- Spatial_Integrated.spe_bysamples[[sample]]
    p <- plotSpatialScatterpie(
        x = Spatial_Integrated.i,
        y = mat.i,
        cell_types = colnames(mat.i),
        img = FALSE,
        scatterpie_alpha = 1,
        pie_scale = 0.4) +
        scale_fill_manual(
        values = pal,
        breaks = names(pal))
    png(paste("Figures/SPOTlight_scatterpie_", refdata, "_", celltype2use, "_", sample, ".png", sep = ""), res = 300, height = 1500, width = 2500)
    print(p + ggplot2::ggtitle(sample))
    dev.off()
}

# Calculate the cell type composition and compare Healthy and IPF samples
# Cell types distribution for each sample (ratio)
CP <- as.matrix(table(Idents(Spatial_Integrated), Spatial_Integrated$group))
CP <- t(t(CP)/colSums(CP)) * 100
CP <- round(CP,1)
CP <- t(CP)
col_fun <-  circlize::colorRamp2(c(0, 100), c("gray99", "red"))
hm <- Heatmap(CP, col = col_fun,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = '#',
        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.0f", CP[i, j]), x, y, gp = gpar(fontsize = 10))})
png(paste0("Figures/CCDistr_ratio_HealthyVSIPF", refdata, "_", celltype2use, ".png"), width = 2000, height = 1000, res = 300)
draw(hm)
dev.off()
write.csv(CP, file = paste0("Results/celltype_composition_ratio_HealthyVSIPF", refdata, "_", celltype2use, ".csv"), quote = F)

# Cell types distribution for each sample (count)
CP <- as.matrix(table(Idents(Spatial_Integrated), Spatial_Integrated$group))
CP <- t(CP)
col_fun <-  circlize::colorRamp2(c(0, mean(CP)), c("gray99", "red"))
hm <- Heatmap(CP, col = col_fun,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        name = '#',
        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.0f", CP[i, j]), x, y, gp = gpar(fontsize = 10))})
png(paste0("Figures/CCDistr_counts_HealthyVSIPF", refdata, "_", celltype2use, ".png"), width = 2500, height = 1000, res = 300)
draw(hm)
dev.off()
write.csv(CP, file = paste0("Results/celltype_composition_count_HealthyVSIPF", refdata, "_", celltype2use, ".csv"), quote = F)

# scProportion Test
prop_test <- sc_utils(Spatial_Integrated)
prop_test <- permutation_test(prop_test, cluster_identity = paste0(refdata, "_", celltype2use, "_SPOTlight"), sample_1 = "Healthy", sample_2 = "IPF", sample_identity = "group")
png(paste0("Figures/scProportionTest_SPOTlight_", refdata, "_", celltype2use, ".png"), res = 300, width = 2000, height = 1250)
print(permutation_plot(prop_test))
dev.off()

print(session_info())


