# DEG, ComplexHeatmap, GO&KEGG by clusterProfiler

# Snapshot and libraries
snapshot <- "/path_to_R_snapshot/"
.libPaths(snapshot)
# BiocManager::install("glmGamPoi", lib=snapshot)
# install.packages("anndata", lib=snapshot)
# reticulate::install_miniconda()
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
# devtools::install_github("https://github.com/MarcElosua/SPOTlight", ref = "bioc_rcpp", lib = snapshot)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(ggcorrplot)
library(devtools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)

gg_color_hue <- function(n) {
        hues <- seq(15, 375, length = n+1)
        hcl(h=hues, l=65, c=100)[1:n]
}

# Set working directory
setwd("/your_working_directory/")

### Load the integrated spatial dataset with cell type annotation
cat("Cell type annotation with SPOTlight!\n")
cat("Loading the integrated spatial dataset and lung reference.\n")
dim_sel <- 20
res_sel <- 0.5
refdata <- "LungCellAtlas"
celltype2use <- "ann_level_3"
Spatial_Integrated <- readRDS(paste0("Results/Spatial_Integrated_", refdata, "_", celltype2use, "_Annotated.rds"))

### Add clinical information
# Sex information
Spatial_Integrated$sex <- ""
Spatial_Integrated$sex[Spatial_Integrated$orig.ident %in% c("Lung_23", "Lung_25", "Lung_62", "Lung_7")] <- "M"
Spatial_Integrated$sex[Spatial_Integrated$orig.ident %in% c("Lung_82", "Lung_57")] <- "F"

# PCA plot of sex information
Spatial_Integrated_4pca <- Spatial_Integrated
Idents(Spatial_Integrated_4pca) <- Spatial_Integrated_4pca$orig.ident
png("Figures/pca_plot_allsamples.png", width = 1500, height = 1250, res = 300)
print(DimPlot(Spatial_Integrated_4pca, reduction = "pca", group.by="orig.ident"))
dev.off()
png("Figures/pca_plot_sex.png", width = 1500, height = 1250, res = 300)
print(DimPlot(Spatial_Integrated_4pca, reduction = "pca", group.by="sex"))
dev.off()
png("Figures/pca_plot_sex_split.png", width = 2000, height = 1250, res = 300)
print(DimPlot(Spatial_Integrated_4pca, reduction = "pca", split.by="sex"))
dev.off()

### DEG identifications
cat("DEG identifications...\n")
# Spatial_Integrated <- readRDS(paste0("Results/Spatial_Integrated_", refdata, "_Annotated.rds"))

# Set default assay as SCT
DefaultAssay(Spatial_Integrated) <- "SCT"

# DEGs for ann lv 3
# Set Idents to ann lv 3
Idents(Spatial_Integrated) <- Spatial_Integrated$LungCellAtlas_ann_level_3_SPOTlight
# DEGs for each cell type when comparing IPF vs Healthy
degs_list <- list() # List to store DEGs for each cell type
Spatial_Integrated <- PrepSCTFindMarkers(object = Spatial_Integrated) # Run PrepSCTFindMarkers before running `FindMarkers()`
for (cell_type in unique(Spatial_Integrated$LungCellAtlas_ann_level_3_SPOTlight)) {
    cat("Processing DEGs for", cell_type, "...\n")
    cell_data <- subset(Spatial_Integrated, subset = LungCellAtlas_ann_level_3_SPOTlight == cell_type)
    if (sum(table(cell_data$group)>=50)<2) {
        next
    }
    if (sum(table(cell_data$group)<3)) {
        cat("Cell groups has fewer than 3 cells for", cell_type, "...\n")
        table(cell_data$group)
        next # Skip if either IPF or Healthy samples have less than 50 cells
    }
    # cell_data <- PrepSCTFindMarkers(object = cell_data) # Run PrepSCTFindMarkers before running `FindMarkers()`
    degs <- FindMarkers(cell_data, ident.1 = "IPF", ident.2 = "Healthy", assay = "SCT", group.by = 'group', min.pct = 0.1, logfc.threshold = 0.2, recorrect_umi = FALSE)   # find the difference genes
    degs <- degs[degs$p_val_adj<0.05, ] # Filter degs by adjusted p value less than 0.05
    if (nrow(degs)>0) {
        degs <- cbind(gene = rownames(degs), cbind(degs, cell_type))
        cat("# DEGs for ", cell_type, ":", nrow(degs), "\n")
        cat("# Up regulated DEGs for ", cell_type, ":", sum(degs$avg_log2FC>0), "\n")
        cat("# Downregulated DEGs for ", cell_type, ":", sum(degs$avg_log2FC<0), "\n")
    } else {
        degs <- cbind(gene=character(nrow(degs)), cbind(degs, cell_type=character(nrow(degs))))
        cat("# DEGs for ", cell_type, ":", nrow(degs), "\n")
    }
    degs_list[[cell_type]] <- degs
}
# Compare deg matrix
degs_all <- do.call(rbind, degs_list)
write.csv(degs_all, file = "Results/DEGs_all.csv", quote = F, row.names = F)

# Create heatmaps
degs_all <- read.csv("Results/DEGs_all.csv")
group.combined_hm <- Spatial_Integrated
Idents(group.combined_hm) <- group.combined_hm$group
DEG_stats <- c()
for (celltype in unique(degs_all$cell_type)) {
    cat("Processing", celltype, "...\n")
    # DEGs for the cell type
    # celltype <- unique(degs_all$cell_type)[i]
    DEGs <- degs_all[degs_all$cell_type == celltype, ]
    # DEGs <- DEGs[abs(DEGs$avg_log2FC)>7, ]
    DEGs <- DEGs[order(DEGs$avg_log2FC, decreasing = T), ]
    if (nrow(DEGs)<5) {
        next # Skip if number of DEGs is less than 5
    }
    DEGs_cate <- rep("", nrow(DEGs))
    names(DEGs_cate) <- DEGs$genes
    DEGs_cate[DEGs$avg_log2FC>0] <- "up"
    DEGs_cate[DEGs$avg_log2FC<0] <- "dn"
    
    # # Heatmap
    # default_colors <- Seurat::PurpleAndYellow()
    # # Make the color scheme symmetrical
    # # max_val <- max(abs(group.combined_hm@assays$RNA@scale.data))
    # # breaks <- seq(-2, 2, length.out = length(default_colors))
    group.combined_hm.i <- group.combined_hm[, group.combined_hm$LungCellAtlas_ann_level_3_SPOTlight == celltype]
    hm2plot <- group.combined_hm.i[DEGs$gene, ]
    VariableFeatures(hm2plot) <- rownames(hm2plot)
    
    # p <- DoHeatmap(hm2plot, size = 4, slot = "scale.data") + theme(axis.text.y = element_text(size = 12)) + scale_fill_gradient2(low = "purple", mid = "black", high = "yellow", midpoint = 0)
    # png(paste0("Figures/Heatmap_DEGs_", celltype, ".png"), res = 300, height = 2250, width = 2000)
    # print(p)
    # dev.off()

    # Extract the data
    df <- as.matrix(hm2plot@assays$SCT@data[DEGs$gene, ])
    png(paste0("Figures/Histogram_DEGs_exprdistr_", celltype, ".png"), res = 300, height = 1500, width = 1750)
    hist(df, 50, col = "blue", border = F)
    dev.off()
    
    # Create the row annotation using DEGs_cate and create the column annotation: 1. IPF vs Healthy, 2. Samples
    row_anno <- DEGs_cate
    row_anno <- factor(row_anno, levels = c("up", "dn"))
    anno_category <- factor(hm2plot$group)
    anno_samples <- factor(hm2plot$orig.ident)
    anno_samples_col <- gg_color_hue(length(unique(anno_samples)))
    names(anno_samples_col) <- unique(anno_samples)
    anno_colors <- list(
        Category = c("IPF" = gg_color_hue(2)[1], "Healthy" = gg_color_hue(2)[2]),
        Samples = anno_samples_col,
        DEGs = c("up" = "yellow", "dn" = "purple"))
    ra <- rowAnnotation(DEGs = row_anno, col = list(DEGs = anno_colors$DEGs))
    ha <- HeatmapAnnotation(Category = anno_category, Samples = anno_samples, col = anno_colors)

    # Plot the heatmap
    val_limit <- mean(df)+2*sd(df)
    p <- Heatmap(df,
        name = "Expression",
        top_annotation = ha,
        right_annotation = ra,
        row_split = row_anno,
        column_split = anno_category,
        col = colorRamp2(c(0, val_limit), c("black", "yellow")),
        show_column_names = F,
        show_row_names = F,
        cluster_rows = F,
        cluster_columns = F)
    png(paste0("Figures/Heatmap_DEGs_", celltype, ".png"), res = 300, height = 1500, width = 1750)
    print(p)
    dev.off()

    # Update DEG_stats.i
    DEG_stats.i <- c(celltype, ncol(group.combined_hm.i), nrow(DEGs), sum(DEGs_cate=="up"), sum(DEGs_cate=="dn"))
    DEG_stats <- rbind(DEG_stats, DEG_stats.i)
}
DEG_stats <- as.data.frame(DEG_stats)
rownames(DEG_stats) <- NULL
colnames(DEG_stats) <- c("celltype", "# cells", "# DEGs", "# DEGs up", "# DEGs dn")
DEG_stats <- DEG_stats[order(as.numeric(DEG_stats[, "# DEGs"]), decreasing = T), ]
DEG_stats

# DEG volcano plot by sex information
degs_all <- read.csv("Results/DEGs_all.csv")
Spatial_Integrated_volcanoplot <- Spatial_Integrated
degs_sex <- FindMarkers(Spatial_Integrated_volcanoplot, ident.1 = "M", ident.2 = "F", assay = "SCT", group.by = "sex", min.pct = 0.01, logfc.threshold = 0.01, recorrect_umi = FALSE)
# genes_sex <- intersect(rownames(degs_sex), degs_all$gene)
# degs_sex <- degs_sex[genes_sex, ]
# Create the plot
p <- ggplot(degs_sex, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) +  # scatter plot
  theme_minimal() +          # clean theme
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 adj P-value"
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 300), xlim = c(-1, 1))
png("Figures/volcano_plot_sex_degs.png", res = 300, height = 1250, width = 1250)
print(p)
dev.off()

# Filter DEGs
degs_sex_filtered <- degs_sex[abs(degs_sex$avg_log2FC)>0.2, ]
degs_sex_filtered <- degs_sex_filtered[degs_sex_filtered$pct.1>0.1 | degs_sex_filtered$pct.2>0.1, ]
degs_sex_filtered <- degs_sex_filtered[degs_sex_filtered$p_val_adj<0.05, ]
length(intersect(rownames(degs_sex_filtered), degs_all$genes))

### Pathway enrichment analysis
cat("Pathway enrichment analysis...\n")
# Source R function
source("/path_to_R_code/functions.R")

# Create result and figure direcotry for enriched pathways
dir.create("Results/Pathway_enrich")
dir.create("Figures/Pathway_enrich")

# Download the KEGG files
code_dir_name <- "/path_to_code_dir/"
prepare_KEGG <- function(species, KEGG_Type="KEGG", keyType="kegg") {
    keggg_db_file <-paste0(code_dir_name,"kegg_data/kegg_",species,"_20230703.RDS") 
    if(file.exists(keggg_db_file))
    {
        kegg <- readRDS(keggg_db_file)
    }else
    {
        kegg <- clusterProfiler::download_KEGG(species, KEGG_Type, keyType)
    }
    build_Anno(kegg$KEGGPATHID2EXTID,
               kegg$KEGGPATHID2NAME)
}

# Default number of pathway results to show on the plot
default_top_pathways <- 20

# Pathway analysis - KEGG
pathway_celltypes <- names(which(unlist(lapply(degs_list, "nrow"))>4))
for (i in 1:length(pathway_celltypes)) {
    # Cell type i and degs
    celltype.i <- pathway_celltypes[i]
    cat("Pathway enrichment analysis for", celltype.i, "...\n")
    symbol.i <- rownames(degs_list[[celltype.i]])
    entrezID.i = bitr(symbol.i, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID

    # Run KEGG pathway analyses
    kegg_result.i <- enrichKEGG_custom(gene = entrezID.i, organism = "human", pvalueCutoff = 0.05)
    kegg_result_df.i <- as.data.frame(kegg_result.i)
    #Convert Gene Id to symbol
    if(nrow(kegg_result_df.i) > 1){
        kegg_result_df.i$geneSymbols <- unlist(lapply(kegg_result_df.i$geneID, function(geneID){
        geneIDs <- unlist(stringr::str_split(geneID,"/"))
            DE_gene_symbol = bitr(geneIDs, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")$SYMBOL
            return(paste(DE_gene_symbol,collapse = "/"))
        }))
    }
    write.csv(kegg_result_df.i, file = paste0("Results/Pathway_enrich/KEGG_enrich_", celltype.i, ".csv"))

    # Write pngs using information from the longest text name result since that will impact the graphic
    # KEGG pathways
    if (length(kegg_result.i$Description)!=0) {
        max_desc <- max(sapply(kegg_result.i$Description,nchar))
        kegg_result.i <- as.data.frame(kegg_result.i)
        kegg_result.i$Description <- gsub(" - .*", "", kegg_result.i$Description)
        kegg_result.i <- new("enrichResult", result = kegg_result.i)
        if(length(kegg_result.i$Description)<default_top_pathways){
            height_calc <- (60+(length(kegg_result.i$Description)*10))
        } else {
            max_desc <- max(sapply(kegg_result.i$Description[1:default_top_pathways],nchar))
            height_calc <- (60+(default_top_pathways)*10)
        }
        if (max_desc <= 50) {
            width_calc <- (300+(max_desc)*5)
        } else {
            width_calc = 550
        }
        if (height_calc > width_calc) {
            width_calc = height_calc
        }
        png(paste0("Figures/Pathway_enrich/KEGG_enrich_", celltype.i, ".png"), width=width_calc*4, height=height_calc*4, res=300)
        print(graphics::barplot(kegg_result.i, x="GeneRatio",
            showCategory=default_top_pathways,title = paste("KEGG","",sep=""),font.size=8) +
            scale_x_continuous(, expand = expansion(mult = c(0, .1))) +
            scale_y_discrete(labels=function(x)stringr::str_trunc(x, 50)) +
            theme(plot.title = element_text(size=8), legend.title = element_text(size=8),
                legend.text = element_text(size=8),
                legend.key.size=unit(((0.25+(min(default_top_pathways,length(kegg_result.i$Description))*0.03) +
                     (max_desc*0.0001))), units = "cm")))
        dev.off()
    }
}

# GO analysis - CC
for (i in 1:length(pathway_celltypes)) {
    # Cell type i and degs
    celltype.i <- pathway_celltypes[i]
    cat("Pathway enrichment analysis for", celltype.i, "...\n")
    symbol.i <- rownames(degs_list[[celltype.i]])
    entrezID.i = bitr(symbol.i, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID

    # Run GO analysis
    go_cc_all.i <- enrichGO(gene = entrezID.i, OrgDb = "org.Hs.eg.db", ont = "CC", pvalueCutoff = 0.05, readable = TRUE)
    write.csv(as.data.frame(go_cc_all.i), file = paste0("Results/Pathway_enrich/GO_CC_enrich_", celltype.i, ".csv"))

    # GO CC
    if (length(go_cc_all.i$Description)!=0) {
        max_desc <- max(sapply(go_cc_all.i$Description,nchar))
        if(length(go_cc_all.i$Description)<default_top_pathways){
            height_calc <- (60+(length(go_cc_all.i$Description)*10))
        } else {
            max_desc <- max(sapply(go_cc_all.i$Description[1:default_top_pathways],nchar))
            height_calc <- (60+(default_top_pathways)*10)
        }
        if (max_desc <= 50) {
            width_calc <- (300+(max_desc)*5)
        } else {
            width_calc = 550
        }
        if (height_calc > width_calc) {
            width_calc = height_calc
        }
        png(paste0("Figures/Pathway_enrich/GO_CC_enrich_", celltype.i, ".png"), width=width_calc*4, height=height_calc*4, res=300)
        print(graphics::barplot(go_cc_all.i, x="GeneRatio",
            showCategory=default_top_pathways,title = paste("GO_CC","",sep=""),font.size=8) +
            scale_x_continuous(, expand = expansion(mult = c(0, .1))) +
            scale_y_discrete(labels=function(x)stringr::str_trunc(x, 50)) +
            theme(plot.title = element_text(size=8), legend.title = element_text(size=8),
                legend.text = element_text(size=8),
                legend.key.size=unit(((0.25+(min(default_top_pathways,length(go_cc_all.i$Description))*0.03) +
                     (max_desc*0.0001))), units = "cm")))
        dev.off()
    }
}

# GO analysis - BP
for (i in 1:length(pathway_celltypes)) {
    # Cell type i and degs
    celltype.i <- pathway_celltypes[i]
    cat("Pathway enrichment analysis for", celltype.i, "...\n")
    symbol.i <- rownames(degs_list[[celltype.i]])
    entrezID.i = bitr(symbol.i, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID

    # Run GO analysis
    go_bp_all.i <- enrichGO(gene = entrezID.i, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
    write.csv(as.data.frame(go_bp_all.i), file = paste0("Results/Pathway_enrich/GO_BP_enrich_", celltype.i, ".csv"))
    
    # GO BP
    if (length(go_bp_all.i$Description)!=0) {
        max_desc <- max(sapply(go_bp_all.i$Description,nchar))
        if(length(go_bp_all.i$Description)<default_top_pathways){
            height_calc <- (60+(length(go_bp_all.i$Description)*10))
        } else {
            max_desc <- max(sapply(go_bp_all.i$Description[1:default_top_pathways],nchar))
            height_calc <- (60+(default_top_pathways)*10)
        }
        if (max_desc <= 50) {
            width_calc <- (300+(max_desc)*5)
        } else {
            width_calc = 550
        }
        if (height_calc > width_calc) {
            width_calc = height_calc
        }
        png(paste0("Figures/Pathway_enrich/GO_BP_enrich_", celltype.i, ".png"), width=width_calc*4, height=height_calc*4, res=300)
        print(graphics::barplot(go_bp_all.i, x="GeneRatio",
            showCategory=default_top_pathways,title = paste("GO_BP","",sep=""),font.size=8) +
            scale_x_continuous(, expand = expansion(mult = c(0, .1))) +
            scale_y_discrete(labels=function(x)stringr::str_trunc(x, 50)) +
            theme(plot.title = element_text(size=8), legend.title = element_text(size=8),
                legend.text = element_text(size=8),
                legend.key.size=unit(((0.25+(min(default_top_pathways,length(go_bp_all.i$Description))*0.03) +
                     (max_desc*0.0001))), units = "cm")))
        dev.off()
    }
}

# GO analysis - MF
for (i in 1:length(pathway_celltypes)) {
    # Cell type i and degs
    celltype.i <- pathway_celltypes[i]
    cat("Pathway enrichment analysis for", celltype.i, "...\n")
    symbol.i <- rownames(degs_list[[celltype.i]])
    entrezID.i = bitr(symbol.i, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID

    # Run GO analysis
    go_mf_all.i <- enrichGO(gene = entrezID.i, OrgDb = "org.Hs.eg.db", ont = "MF", pvalueCutoff = 0.05, readable = TRUE)
    write.csv(as.data.frame(go_mf_all.i), file = paste0("Results/Pathway_enrich/GO_MF_enrich_", celltype.i, ".csv"))

    # GO MF
    if (length(go_mf_all.i$Description)!=0) {
        max_desc <- max(sapply(go_mf_all.i$Description,nchar))
        if(length(go_mf_all.i$Description)<default_top_pathways){
            height_calc <- (60+(length(go_mf_all.i$Description)*10))
        } else {
            max_desc <- max(sapply(go_mf_all.i$Description[1:default_top_pathways],nchar))
            height_calc <- (60+(default_top_pathways)*10)
        }
        if (max_desc <= 50) {
            width_calc <- (300+(max_desc)*5)
        } else {
            width_calc = 550
        }
        if (height_calc > width_calc) {
            width_calc = height_calc
        }
        png(paste0("Figures/Pathway_enrich/GO_MF_enrich_", celltype.i, ".png"), width=width_calc*4, height=height_calc*4, res=300)
        print(graphics::barplot(go_mf_all.i, x="GeneRatio",
            showCategory=default_top_pathways,title = paste("GO_MF","",sep=""),font.size=8) +
            scale_x_continuous(, expand = expansion(mult = c(0, .1))) +
            scale_y_discrete(labels=function(x)stringr::str_trunc(x, 50)) +
            theme(plot.title = element_text(size=8), legend.title = element_text(size=8),
                legend.text = element_text(size=8),
                legend.key.size=unit(((0.25+(min(default_top_pathways,length(go_mf_all.i$Description))*0.03) +
                     (max_desc*0.0001))), units = "cm")))
        dev.off()
    }
}

# Print out the R session details
print(session_info())




