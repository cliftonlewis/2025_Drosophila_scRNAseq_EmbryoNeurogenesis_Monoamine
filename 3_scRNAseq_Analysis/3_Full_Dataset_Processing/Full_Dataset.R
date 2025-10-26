# ============================
# Integrated Analysis Script
# ============================

# ────────────────────────────────────────────────────
# 1. Environment setup
# ────────────────────────────────────────────────────
.libPaths('~/R/x86_64-pc-linux-gnu-library/4.5/')
options(future.globals.maxSize = 40000 * 1024^2)
setwd("/DataDrives/Drive2/Clifton/R_Projects/2025_Drosophila_scRNAseq_MonoamineSpecification/ANALYSIS/Step_2_Initial_preprocessing/Step_2.3_Full_Dataset_preprocessing/")

# ────────────────────────────────────────────────────
# 2. Load required packages
# ────────────────────────────────────────────────────
# --- scRNA-seq analysis ---
library(scCustomize)
library(Seurat)
library(SeuratObject)
library(sctransform)
library(clustree)
library(glmGamPoi)
library(DoubletFinder)
library(harmony)
library(monocle3)
library(SeuratWrappers)
library(SeuratExtend)
library(org.Dm.eg.db)
library(biomaRt)
library(RAPToR)
library(drosoRef)

# --- Data science & plotting ---
library(tidyverse)
library(readr)
library(stringr)
library(dplyr)
library(qs)
library(magrittr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dittoSeq)

# --- Annotation database ---
organism <- "org.Dm.eg.db"
library(organism, character.only = TRUE)

# ────────────────────────────────────────────────────
# 3. Source custom functions
# ────────────────────────────────────────────────────
file_path2 <- paste0("/DataDrives/Drive2/Clifton/R_Projects/2025_Drosophila_scRNAseq_MonoamineSpecification/ANALYSIS/Step_2_Initial_preprocessing")

source(file.path(file_path2, "Functions/sccustomize_qc.R"))
source(file.path(file_path2, "Functions/create_seurat_object.R"))
source(file.path(file_path2, "Functions/add_metadata.R"))
source(file.path(file_path2, "Functions/perform_qc_and_plots.R"))
source(file.path(file_path2, "Functions/perform_qc_filtering.R"))
source(file.path(file_path2, "Functions/run_normalisation.R"))
source(file.path(file_path2, "Functions/cell_cycle_scoring.R"))
source(file.path(file_path2, "Functions/generate_var_features_plot.R"))
source(file.path(file_path2, "Functions/run_pca.R"))
source(file.path(file_path2, "Functions/determine_pcs.R"))
source(file.path(file_path2, "Functions/run_dimred.R"))
source(file.path(file_path2, "Functions/perform_clustree.R"))
source(file.path(file_path2, "Functions/perform_clustree_harmony.R"))
source(file.path(file_path2, "Functions/doubletfinder.R"))
source(file.path(file_path2, "Functions/fishers_test.R"))  

# ────────────────────────────────────────────────────
# 4. Load data & define global parameters
# ────────────────────────────────────────────────────
# Load datasets
hrs_00_03 <- qread("../Step_2.2_Timepoint_preprocessing/Data/Timepoint_Datasets/hrs_00_03/3_hrs_00_03_harmonyintegratedannotated.qs")
hrs_02_05 <- qread("../Step_2.2_Timepoint_preprocessing/Data/Timepoint_Datasets/hrs_02_05/3_hrs_02_05_harmonyintegratedannotated.qs")
hrs_04_07 <- qread("../Step_2.2_Timepoint_preprocessing/Data/Timepoint_Datasets/hrs_04_07/3_hrs_04_07_harmonyintegratedannotated.qs")
hrs_05_08 <- qread("../Step_2.2_Timepoint_preprocessing/Data/Timepoint_Datasets/hrs_05_08/3_hrs_05_08_harmonyintegratedannotated.qs")
hrs_06_10 <- qread("../Step_2.2_Timepoint_preprocessing/Data/Timepoint_Datasets/hrs_06_10/3_hrs_06_10_harmonyintegratedannotated.qs")
hrs_09_13 <- qread("../Step_2.2_Timepoint_preprocessing/Data/Timepoint_Datasets/hrs_09_13/3_hrs_09_13_harmonyintegratedannotated.qs")
hrs_12_17 <- qread("../Step_2.2_Timepoint_preprocessing/Data/Timepoint_Datasets/hrs_12_17/3_hrs_12_17_harmonyintegratedannotated.qs")
hrs_16_22 <- qread("../Step_2.2_Timepoint_preprocessing/Data/Timepoint_Datasets/hrs_16_22/3_hrs_16_22_harmonyintegratedannotated.qs")

# Merge all Seurat objects
merged_seurat <- merge(
  hrs_00_03, 
  y = list(hrs_02_05, hrs_04_07, hrs_05_08, hrs_06_10, hrs_09_13, hrs_12_17, hrs_16_22)
)
rm(hrs_00_03)
rm(hrs_02_05)
rm(hrs_04_07)
rm(hrs_05_08)
rm(hrs_06_10)
rm(hrs_09_13)
rm(hrs_12_17)
rm(hrs_16_22)

# Check if 'celltype' exists in the metadata
if (!"celltype" %in% colnames(merged_seurat@meta.data)) {
  stop("The column 'celltype' does not exist in the merged dataset metadata.")
}

# Replace blank or NA values in 'celltype' with 'Unknown'
merged_seurat@meta.data$celltype[is.na(merged_seurat@meta.data$celltype) | 
                                   merged_seurat@meta.data$celltype == ""] <- "Unknown"

# Verify changes
table(merged_seurat@meta.data$celltype)

# Optional: Save the updated Seurat object
qsave(merged_seurat, file = "merged_timepoints.qs")
merged_seurat <- qread(file = "merged_timepoints.qs")
sample_name <- "full_dataset"
pcs <- 100
merged_seurat <- run_normalisation(merged_seurat)
merged_seurat <- JoinLayers(merged_seurat)
merged_seurat <- run_normalisation(merged_seurat)
merged_seurat <- cell_cycle_scoring(sample_object = merged_seurat)
merged_seurat <- generate_variable_features_plot(merged_seurat, sample_name)
# Run PCA
merged_seurat <- RunPCA(object = merged_seurat, features = VariableFeatures(object = merged_seurat), npcs = 200)
# Plot the elbow plot
plot1 <- ElbowPlot(object = merged_seurat, ndims=200)
print(plot1)
ggsave(paste0("Plots/Timepoint_Datasets/", sample_name, "/4_", sample_name,  "_ElbowPlot.png"), plot1, width = 10, height = 6)

# Calculate the percentage of variation associated with each PC
pct <- merged_seurat[["pca"]]@stdev / sum(merged_seurat[["pca"]]@stdev) * 100

# Calculate the cumulative percentage of variation for each PC
cumu <- cumsum(pct)

# Method 1a: Determine the first PC where the cumulative variance exceeds 90% 
# and the current PC contributes less than 5%
co1 <- which(cumu > 90 & pct < 5)[1]

# Method 1b: Determine the first PC where the percent change in variation between 
# consecutive PCs is less than 0.1%
co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1

# Use the minimum of the two calculations to determine the optimal number of PCs
pcs_variance <- min(co1, co2)

# Create a dataframe for plotting the elbow plot
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

# Elbow plot to visualize the percentage of variation explained by each PC
plot1 <- ggplot(plot_df, aes(x = rank, y = pct)) + 
  geom_point() +
  geom_line() +
  geom_vline(xintercept = pcs_variance, linetype = "dashed", color = "red", size = 1) +  # Red dashed line for optimal PC
  geom_text(aes(label = rank), nudge_y = 0.5, size = 3) +
  ggtitle("Elbow Plot: Variance Explained by Principal Components") +
  labs(x = "Principal Components", y = "Percentage of Variance Explained (%)") +
  annotate("text", x = max(plot_df$rank) - 5, y = max(pct) - 10, label = paste("PCs chosen: ", pcs_variance), 
           color = "blue", size = 5, hjust = 1, vjust = 1)  # Label for the chosen PCs in the top right corner
print(plot1)
# Save the elbow plot
ggsave(paste0("Plots/Timepoint_Datasets/", sample_name, "/5A_", sample_name, "_ElbowPlot_significant.png"), plot1, width = 10, height = 6)

print(paste("Optimal number of PCs based on variance explained:", pcs_variance))
# 
# # Run the JackStraw procedure for determining statistically significant PCs
# merged_seurat <- JackStraw(merged_seurat, dims = 200)
# merged_seurat <- ScoreJackStraw(merged_seurat, dims = 1:200)
# 
# # JackStraw plot to visualize the statistical significance of each PC
# pcs_jackstraw <- length(which(merged_seurat[["pca"]]@jackstraw@overall.p.values[, "Score"] < 0.05))
# 
# plot2 <- JackStrawPlot(merged_seurat, dims = 1:200) +
#   annotate("text", x = pcs_jackstraw + 3, y = max(merged_seurat[["pca"]]@jackstraw@overall.p.values[, "Score"]), 
#            label = paste("PCs chosen: ", pcs_jackstraw), color = "blue", size = 5) +  # Annotate chosen PCs
#   ggtitle("JackStraw Plot: Significance of Principal Components") +
#   labs(x = "Principal Components", y = "p-value")
# 
# # Save the JackStraw plot
# ggsave(paste0("Plots/Timepoint_Datasets/", sample_name, "/5B_", sample_name, "_JackStrawPlot_significant.png"), plot2, width = 20, height = 6)
# 
# print(paste("Optimal number of PCs based on JackStraw plot:", pcs_jackstraw))
# 
# # Find neighbors
# merged_seurat <- FindNeighbors(merged_seurat, dims = 1:pcs)
# 
# # Perform clustering
# resolution.range <- seq(from = 0, to = 4, by = 0.1)
# merged_seurat <- FindClusters(merged_seurat,
#                               reduction = "pca", dims = 1:pcs,
#                               resolution = resolution.range,
#                               save.SNN = TRUE,
#                               print.output = FALSE)
# 
# # Generate clustree plot
# plot1 <- clustree(merged_seurat) + theme_light()
# ggsave(paste0("Plots/Timepoint_Datasets/", sample_name, "/6_", sample_name, "_Clustree.png"), plot1, width = 10, height = 15)
# 
# # Generate clustree plot with SC3 stability color-coding
# plot2 <- clustree(merged_seurat, node_colour = "sc3_stability")
# ggsave(paste0("Plots/Timepoint_Datasets/", sample_name, "/7_", sample_name, "_Clustree_sc3_stability.png"), plot2, width = 10, height = 15)

merged_seurat <- FindNeighbors(merged_seurat, dims = 1:pcs)
qsave(merged_seurat, file = "1.qs")
merged_seurat <- qread("1.qs")
merged_seurat <- FindClusters(merged_seurat, reduction = "pca", dims = 1:pcs, resolution = 2.2)
merged_seurat <- run_dimred(merged_seurat, sample_name, pcs)
qsave(merged_seurat, file = paste0("Data/Timepoint_Datasets/", sample_name, "/2_", sample_name, "_preprocessed.qs"))

# Verify the merged object
print(merged_seurat)
head(merged_seurat@meta.data)

full_dataset <- qread(paste0("Data/Timepoint_Datasets/", sample_name, "/2_", sample_name, "_preprocessed.qs"))
# full_dataset <- merged_seurat 

metadata <- full_dataset@meta.data

# Function to get the most frequent celltype for each seurat_cluster
get_most_frequent_label <- function(cluster_id) {
  # Subset the metadata for the specific cluster
  cluster_data <- metadata[metadata$seurat_clusters == cluster_id, ]
  
  # Filter out NA or NULL values in celltype
  cluster_data <- cluster_data[!is.na(cluster_data$celltype) & cluster_data$celltype != "NULL", ]
  
  # Check if we have any valid values left after filtering
  if (nrow(cluster_data) == 0) {
    return("Unknown")  # Assign "Unknown" if no valid cell type exists for the cluster
  }
  
  # Find the most frequent celltype within this cluster
  most_frequent <- names(sort(table(cluster_data$celltype), decreasing = TRUE))[1]
  
  return(most_frequent)
}

# Apply the function to each unique seurat_cluster
most_frequent_labels <- sapply(unique(metadata$seurat_clusters), get_most_frequent_label)

# Create a mapping of seurat_cluster -> most frequent celltype
cluster_to_celltype <- setNames(most_frequent_labels, unique(metadata$seurat_clusters))


# Now create a new column combining seurat_cluster and most frequent celltype
metadata$celltype_merge <- paste0(cluster_to_celltype[as.character(metadata$seurat_clusters)]
)

# Now create a new column combining seurat_cluster and most frequent celltype
metadata$celltype_merge_unique <- paste0(
  metadata$seurat_clusters, "_", cluster_to_celltype[as.character(metadata$seurat_clusters)]
)

# **Corrected Assignment:** Add the updated metadata back to the Seurat object
full_dataset@meta.data$celltype_merge <- metadata$celltype_merge
full_dataset@meta.data$celltype_merge_unique <- metadata$celltype_merge_unique

# Optionally, view the result
table(full_dataset@meta.data$celltype_merge)
table(full_dataset@meta.data$celltype_merge_unique)

# Plot 1: Full Dataset - Annotated by BDGP Putative Cell Type
Idents(full_dataset) <- full_dataset$celltype_merge

# Load your cluster-to-color mapping from a CSV file
# The CSV file should have two columns: "Label" and "Color"
annotation_table <- read.csv("../Supplementary_Data/label2colour_mapping.csv")

# Ensure the Seurat object has cluster identities assigned
clusters_in_seurat <- levels(Idents(full_dataset))  # Get unique clusters in the Seurat object

# Subset the annotation table to include only clusters present in the Seurat object
filtered_annotation_table <- annotation_table[annotation_table$Label %in% clusters_in_seurat, ]

# Create a named vector for colors based on the filtered annotation table
color_mapping_bdgp <- setNames(filtered_annotation_table$Color, filtered_annotation_table$Label)

# Map colors to each cell in the Seurat object
cell_colors <- color_mapping_bdgp[as.character(full_dataset@meta.data$celltype_merge)]

# Add the colors as metadata in the Seurat object
full_dataset@meta.data$celltype_merge_color <- cell_colors

# Verify the new metadata column
head(full_dataset@meta.data)

plot1 <- dittoDimPlot(object = full_dataset, 'celltype_merge', reduction.use = "umap", color.panel = color_mapping_bdgp) +
  theme_minimal() +
  labs(title = "UMAP Plot", subtitle = "Full Dataset Annotated by BDGP Putative Cell Type") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "right")
print(plot1)
ggsave(filename = paste0("Plots/", "Timepoint_Datasets/full_dataset/", "1_full_dataset_annotated_UMAP2D.pdf"), plot = plot1, width = 13, height = 6)

# Plot 2: Full Dataset - Annotated by BDGP Putative Cell Type (Unique)
full_dataset <- SetIdent(full_dataset, value = "celltype_merge_unique")
Idents(full_dataset)  

# Extract the existing color mapping from metadata
color_mapping_bdgp <- setNames(
  unique(full_dataset@meta.data$celltype_merge_color),
  unique(full_dataset@meta.data$celltype_merge)
)

# Create a mapping for `celltype_merge_unique`
# Extract the `celltype_merge` part of `celltype_merge_unique`
unique_labels <- unique(full_dataset@meta.data$celltype_merge_unique)
celltype_parts <- gsub(".*_", "", unique_labels)  # Extract the celltype_bdgp part
color_mapping_unique <- setNames(color_mapping_bdgp[celltype_parts], unique_labels)

plot2 <- dittoDimPlot(object = full_dataset, 'celltype_merge_unique', reduction.use = "umap", do.label=FALSE, color.panel = color_mapping_unique) +
  theme_minimal() +
  labs(title = "UMAP Plot", subtitle = "Full Dataset Annotated by BDGP Putative Cell Type") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "right")
print(plot2)
ggsave(filename = paste0("Plots/Timepoint_Datasets/full_dataset/",  "/2_full_dataset_annotated_unique_UMAP2D.pdf"), plot = plot2, width = 18, height = 6)

# Plot 3: Full Dataset - Annotated by Timepoint
Idents(full_dataset) <- full_dataset$timepoint

# Load your cluster-to-color mapping from a CSV file
# The CSV file should have two columns: "Label" and "Color"
annotation_table <- read.csv("../Supplementary_Data/timepoint2colour_mapping.csv")

# Ensure the Seurat object has cluster identities assigned
clusters_in_seurat <- levels(Idents(full_dataset))  # Get unique clusters in the Seurat object

# Subset the annotation table to include only clusters present in the Seurat object
filtered_annotation_table <- annotation_table[annotation_table$Label %in% clusters_in_seurat, ]

# Create a named vector for colors based on the filtered annotation table
color_mapping_timepoint <- setNames(filtered_annotation_table$Color, filtered_annotation_table$Label)

# Map colors to each cell in the Seurat object
cell_colors <- color_mapping_timepoint[as.character(full_dataset@meta.data$timepoint)]

# Add the colors as metadata in the Seurat object
full_dataset@meta.data$timepoint_color <- cell_colors

# Verify the new metadata column
head(full_dataset@meta.data)

plot3 <- dittoDimPlot(object = full_dataset, 'timepoint', reduction.use = "umap", color.panel = color_mapping_timepoint) +
  theme_minimal() +
  labs(title = "UMAP Plot", subtitle = "Full Dataset Annotated by Timepoint") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "right")
print(plot3)
ggsave(filename = paste0("Plots/Timepoint_Datasets/full_dataset/",   "/3_full_dataset_timepoint_UMAP2D.pdf"), plot = plot3, width = 10, height = 6)

# Plot 4: Full Dataset - Annotated by Cell Cycle Phase
Idents(full_dataset) <- full_dataset$Phase
plot4 <- dittoDimPlot(object = full_dataset, 'Phase', reduction.use = "umap", do.label = TRUE) +
  theme_minimal() +
  labs(title = "UMAP Plot", subtitle = "Full Dataset Annotated by Cell Cycle Phase") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "right")
ggsave(filename = paste0("Plots/Timepoint_Datasets/full_dataset/",   "/4_full_dataset_phase_UMAP2D.pdf"), plot = plot4, width = 10, height = 6)

qsave(full_dataset, file = "Data/Timepoint_Datasets/full_dataset/2A_full_dataset_annotated_w_unique.qs")

library(reticulate)
use_condaenv('sccustomize')
as.anndata(full_dataset, file_path = "Data/Timepoint_Datasets/full_dataset/", "2A_full_dataset_annotated_w_unique.h5ad")

## Neuronal lineage
full_dataset <- qread("Data/Timepoint_Datasets/full_dataset/2A_full_dataset_annotated_w_unique.qs")

# Define your grouped features as a named list
grouped_features <- list(
  "Neuroblasts" = c("trbl", "dpn", "wor", "ase", "mira", "insc", "stg", "tll"),
  "INPs" = c("erm", "ham", "brat"),
  "GMCs" = c("tap", "dap", "pros", "numb"),
  "New-born neurons" = c("Hey", "E(spl)m6-BFM"),
  "Immature neurons" = c("elav", "CadN", "fne"),
  "Mature neurons" = c("brp", "nSyb", "Syn", "Syt1", "Syt4"),
  "Glia" = c("repo", "alrm", "wrapper")
)

plot5 <- DotPlot2(seu = full_dataset, features = grouped_features, group.by = "seurat_clusters", color_scheme = "BuRd")
ggsave(filename = paste0("Plots/Timepoint_Datasets/full_dataset/",   "/5_full_dataset_neuro_markers_seurat_clusters.pdf"), plot = plot5, width = 20, height = 6)

plot5 <- DotPlot2(seu = full_dataset, features = grouped_features, group.by = "celltype_merge_unique", color_scheme = "BuRd")
ggsave(filename = paste0("Plots/Timepoint_Datasets/full_dataset/",   "/5_full_dataset_neuro_markers_celltype_merge_unique.pdf"), plot = plot5, width = 20, height = 8)

plot5 <- DotPlot2(seu = full_dataset, features = grouped_features, group.by = "celltype_merge", color_scheme = "BuRd")
ggsave(filename = paste0("Plots/Timepoint_Datasets/full_dataset/",   "/5_full_dataset_neuro_markers_celltype_merge.pdf"), plot = plot5, width = 10, height = 8)

plot5 <- DotPlot2(seu = full_dataset, features = grouped_features, group.by = "timepoint", color_scheme = "BuRd")
ggsave(filename = paste0("Plots/Timepoint_Datasets/full_dataset/",   "/5_full_dataset_neuro_markers_timepoint.pdf"), plot = plot5, width = 8, height = 6)

# New neuronal annotations

# Ensure the active identities are set to the clusters
Idents(full_dataset) <- full_dataset$seurat_clusters

# (Optional) Check the current levels to verify the cluster names
levels(full_dataset)

# Define the mapping vector for new cluster annotations
new.cluster.ids <- c(
  "0" = "Immature_neurons/Mature_neurons",
  "1" = "Neuroblasts/INPs/GMCs/Immature_neurons",
  "2" = "Neuroblasts/INPs/GMCs/New-born_neurons/Immature_neurons",
  "3" = "Other",
  "4" = "INPs/Immature_neurons/Mature_neurons",
  "5" = "Neuroblasts/GMCs/Immature_neurons",
  "6" = "Neuroblasts/GMCs/Immature_neurons",
  "7" = "Neuroblasts/GMCs/Immature_neurons",
  "8" = "Other",
  "9" = "Neuroblasts/GMCs/Immature_neurons",
  "10" = "Other",
  "11" = "INPs",
  "12" = "Neuroblasts/GMCs/Immature_neurons",
  "13" = "Neuroblasts/INPs/GMCs/Immature_neurons",
  "14" = "Neuroblasts/GMCs",
  "15" = "Other",
  "16" = "Other",
  "17" = "GMCs/Immature_neurons",
  "18" = "INPs",
  "19" = "Other",
  "20" = "Neuroblasts",
  "21" = "Neuroblasts/GMCs",
  "22" = "Neuroblasts/GMCs",
  "23" = "Neuroblasts/GMCs",
  "24" = "Other",
  "25" = "Neuroblasts/GMCs",
  "26" = "Other",
  "27" = "Neuroblasts/GMCs/Immature_neurons",
  "28" = "Neuroblasts/GMCs/Immature_neurons",
  "29" = "Neuroblasts/GMCs/Immature_neurons",
  "30" = "GMCs/Immature_neurons",
  "31" = "Other",
  "32" = "Other",
  "33" = "Other",
  "34" = "Neuroblasts/Immature_neurons",
  "35" = "Other",
  "36" = "Neuroblasts/INPs",
  "37" = "Neuroblasts/GMCs",
  "38" = "GMCs/Glia",
  "39" = "Neuroblasts/INPs/GMCs/New-born_neurons/Immature_neurons",
  "40" = "INPs",
  "41" = "Immature_neurons",
  "42" = "Neuroblasts/GMCs",
  "43" = "Other",
  "44" = "Neuroblasts/GMCs/Immature_neurons",
  "45" = "Neuroblasts/GMCs",
  "46" = "Other",
  "47" = "Other",
  "48" = "GMCs",
  "49" = "Other",
  "50" = "GMCs",
  "51" = "Neuroblasts/INPs/GMCs",
  "52" = "Other",
  "53" = "Other",
  "54" = "INPs",
  "55" = "GMCs",
  "56" = "Immature_neurons/Mature_neurons",
  "57" = "Other",
  "58" = "Glia",
  "59" = "Other",
  "60" = "Neuroblasts",
  "61" = "Other",
  "62" = "Other",
  "63" = "Other",
  "64" = "Other",
  "65" = "Neuroblasts/GMCs/Immature_neurons",
  "66" = "Immature_neurons",
  "67" = "INPs/Immature_neurons/Mature_neurons",
  "68" = "INPs/GMCs/New-born_neurons/Immature_neurons",
  "69" = "Other",
  "70" = "Other",
  "71" = "Immature_neurons/Mature_neurons",
  "72" = "Neuroblasts/GMCs"
)



full_dataset@meta.data$neuronal_annotation_broad <- full_dataset@meta.data$seurat_clusters
Idents(full_dataset) <- full_dataset$neuronal_annotation_broad
levels(full_dataset)

# Rename the active identities in the Seurat object
full_dataset <- RenameIdents(full_dataset, new.cluster.ids)

# Save the new identities into a metadata column for future reference.
full_dataset@meta.data$neuronal_annotation_broad <- Idents(full_dataset)

# Verify the changes by checking the metadata
head(full_dataset@meta.data)
neuronal_cells <- colnames(full_dataset)[full_dataset$neuronal_annotation_broad == c("Neuroblasts/GMCs/Immature_neurons", "Mature neurons", "New-born neurons")]
plot6 <- DimPlot2(full_dataset, cells.highlight = neuronal_cells, cols.highlight = "firebrick1")+
  labs(title = "UMAP Plot", subtitle = "Full Dataset Annotated by Neuronal Annotation") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "right")
print(plot6)
ggsave(filename = paste0("Plots/Timepoint_Datasets/full_dataset/",   "/6_full_dataset_neuronal_annotation_UMAP2D.pdf"), plot = plot6, width = 10, height = 6)
colour_pal <-c("springgreen", "red", "bisque",  "grey",  "darkgoldenrod", "maroon1", "salmon", "purple", "royalblue3", "green4", "brown4", "blue", "dodgerblue", "darkorange", "darkorchid1", "cyan3", "lightpink1", "lightblue1", "lightgreen1", "lightyellow1")
plot7 <- DimPlot2(full_dataset, features = c("neuronal_annotation_broad"), cols = list("neuronal_annotation_broad" = colour_pal))+
  labs(title = "UMAP Plot") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "right")
print(plot7)
ggsave(filename = paste0("Plots/Timepoint_Datasets/full_dataset/",   "/7_full_dataset_neuronal_annotation_UMAP2D.pdf"), plot = plot7, width = 10, height = 6)

qsave(full_dataset, file = "Data/Timepoint_Datasets/full_dataset/2B_full_dataset_annotated_w_unique.qs")

full_dataset <- qread("Data/Timepoint_Datasets/full_dataset/2B_full_dataset_annotated_w_unique.qs")

library(reticulate)
use_condaenv('sccustomize')
as.anndata(full_dataset, file_path = "Data/Timepoint_Datasets/full_dataset/", "2B_full_dataset_annotated_w_unique.h5ad")

meta <- full_dataset@meta.data

meta$neuronal_annotation_broad_unique <- paste(meta$seurat_clusters,
                            meta$neuronal_annotation_broad,
                            sep = "_")

full_dataset <- AddMetaData(full_dataset, meta$neuronal_annotation_broad_unique, col.name = "neuronal_annotation_broad_unique")
Idents(full_dataset) <- "neuronal_annotation_broad_unique"  

table(full_dataset@meta.data$neuronal_annotation_broad_unique)

# Original neuronal subset
mature_neuronal_subset <- subset(full_dataset, neuronal_annotation_broad %in% c("Neuroblasts/INPs/GMCs/Immature_neurons", "Immature_neurons/Mature_neurons", "INPs/Immature_neurons/Mature_neurons", "GMCs/Glia", "Neuroblasts/INPs/GMCs"))

# # Interactive dim plot with hover info
# p <- dittoDimPlot(
#   full_dataset,
#   var = "neuronal_annotation_broad_unique",
#   do.hover = TRUE,
#   hover.data = c("neuronal_annotation_broad_unique", "seurat_clusters") 
# )
# p <- DimPlot(full_dataset, reduction = "umap", group.by = "neuronal_annotation_broad_unique", label = TRUE)
# picked <- CellSelector(p)  # draw a polygon around the cluster in the plot window
# length(picked); head(picked)

# p

# # Create a subset 
# mature_neuronal_subset <- subset(full_dataset, neuronal_annotation_broad %in% c("Immature_neurons/Mature_neurons", "Neuroblasts/INPs/GMCs/Immature_neurons", "Neuroblasts/INPs/GMCs/New-born_neurons/Immature_neurons", "INPs/Immature_neurons/Mature_neurons", "Neuroblasts/GMCs/Immature_neurons", "INPs", "Neuroblasts/GMCs", "GMCs/Immature_neurons", "Neuroblasts", "Neuroblasts/Immature_neurons", "Neuroblasts/INPs", "Immature_neurons", "GMCs", "Neuroblasts/INPs/GMCs", "INPs/GMCs/New-born_neurons/Immature_neurons"))
# DimPlot(object = mature_neuronal_subset, group.by = "neuronal_annotation_broad")

# # mature_neuronal_subset <- subset(full_dataset, neuronal_annotation_broad %in% c("Neuroblasts/INPs/GMCs/New−born_neurons/Immature_neurons", "Neuroblasts/INPs/GMCs/Immature_neurons", "Neuroblasts/INPs/GMCs", "GMCs/Glia", "INPs/Immature_neurons/Mature_neurons", "Immature_neurons/Mature_neurons"))

# Verify the unique neuronal annotation labels in the subset
table(mature_neuronal_subset@meta.data$neuronal_annotation_broad)

# Plot the UMAP for the neuronal subset using DimPlot2 from seuratextend.
library(RColorBrewer)
# Number of categories in your data
n <- length(unique(mature_neuronal_subset$neuronal_annotation_broad))
# Build a combined palette from actual colors
combined_palette <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(12, "Set3")
)
# Take only as many colors as needed
combined_palette <- combined_palette[1:n]
# Pass this to your plot
plot8 <- DimPlot2(
  mature_neuronal_subset,
  features = "neuronal_annotation_broad",
  cols = combined_palette
) +
  labs(title = "UMAP Plot") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    legend.position = "right"
  )
print(plot8)
ggsave(filename = paste0("Plots/Timepoint_Datasets/mature_neuronal_subset/",   "/1_mature_neuronal_subset_neuronal_annotation_UMAP2D.pdf"), plot = plot8, width = 10, height = 6)

qsave(mature_neuronal_subset, file = "Data/Timepoint_Datasets/mature_neuronal_subset/1_mature_neuronal_subset_preannotated.qs")








# ===========================================================================================================================================
## Gene ontology

# STEP 1: Load the full dataset
full_dataset <- qread("Data/Timepoint_Datasets/full_dataset/2A_full_dataset_annotated_w_unique.qs")

# STEP 2: Get all genes from assays
assays_to_annotate <- c("spliced", "unspliced", "ambiguous", "RNA")
all_genes <- unique(unlist(lapply(assays_to_annotate, function(x) rownames(full_dataset[[x]]))))

# STEP 3: Query biomaRt for gene annotations
ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
annotations <- getBM(
  attributes = c("external_gene_name", "entrezgene_id", "flybase_gene_id"),
  filters = "external_gene_name",
  values = all_genes,
  mart = ensembl
)

# STEP 4: Make a lookup table for annotations
annotations_df <- annotations %>%
  distinct(external_gene_name, .keep_all = TRUE) %>%
  column_to_rownames("external_gene_name")

# STEP 5: Identify top markers using Seurat's FindAllMarkers
top_markers <- FindAllMarkers(
  object = full_dataset,
  only.pos = TRUE,
  group.by = "celltype_merge"
)

# Check if the 'gene' column already exists;
# if not, add it from rownames.
if (!("gene" %in% colnames(top_markers))) {
  top_markers <- top_markers %>% rownames_to_column("gene")
}

top_markers_renamed <- top_markers %>%
  left_join(annotations, by = c("gene" = "external_gene_name")) %>%
  rename(Entrez_ID = entrezgene_id)

# Filter markers by avg_log2FC and adjusted p-value
top_markers.sig <- top_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.25, p_val_adj <= 0.01) %>%
  ungroup()

# STEP 6: Perform GO enrichment for each cluster and collect results
all_genes <- unique(top_markers.sig$gene)
go_results <- list()

for (cluster_id in unique(top_markers.sig$cluster)) {
  cluster_specific_genes <- top_markers.sig %>%
    dplyr::filter(cluster == cluster_id) %>%
    dplyr::pull(gene)
  
  ego <- clusterProfiler::enrichGO(gene = cluster_specific_genes,
                                   universe = all_genes,
                                   OrgDb = org.Dm.eg.db,
                                   keyType = 'SYMBOL',
                                   ont = "BP", # You can also run for "MF" and "CC" separately
                                   pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)
  
  go_results[[paste0("Cluster_", cluster_id)]] <- ego
  print(paste("GO analysis done for Cluster", cluster_id))
}


p2 <- clusterProfiler::compareCluster(
  geneCluster = gene ~ cluster,
  data = top_markers.sig,
  fun = "enrichGO",
  universe = unique(top_markers.sig$gene),
  OrgDb = org.Dm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",             # Change to "MF" or "CC" if needed
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

clusterProfiler::dotplot(p2, showCategory = 1)


# STEP 7: Function to Calculate Z-scores for GO Terms
calculate_go_zscores_from_enrichment <- function(go_results) {
  all_go_terms <- unique(unlist(lapply(go_results, function(x) if (!is.null(x) && nrow(x@result) > 0) x@result$Description else NULL)))
  cluster_ids <- names(go_results)
  zscore_matrix <- matrix(NA, nrow = length(all_go_terms), ncol = length(cluster_ids))
  rownames(zscore_matrix) <- all_go_terms
  colnames(zscore_matrix) <- cluster_ids
  
  for (i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[i]
    ego_result <- go_results[[cluster_id]]
    if (!is.null(ego_result) && nrow(ego_result@result) > 0) {
      # Calculate Z-score based on the enrichment score (if available)
      # A simple Z-score can be calculated relative to the mean and SD of enrichment scores
      enrichment_scores <- ego_result@result$enrichmentFactor
      mean_enrichment <- mean(enrichment_scores, na.rm = TRUE)
      sd_enrichment <- sd(enrichment_scores, na.rm = TRUE)
      
      for (go_term in all_go_terms) {
        term_row <- ego_result@result[ego_result@result$Description == go_term, ]
        if (nrow(term_row) > 0) {
          if (sd_enrichment > 0) {
            zscore <- (term_row$enrichmentFactor - mean_enrichment) / sd_enrichment
            zscore_matrix[go_term, cluster_id] <- zscore
          } else {
            zscore_matrix[go_term, cluster_id] <- 0 # Avoid division by zero if all enrichment factors are the same
          }
        } else {
          zscore_matrix[go_term, cluster_id] <- NA # GO term not enriched in this cluster
        }
      }
    } else {
      zscore_matrix[, cluster_id] <- NA # No GO results for this cluster
    }
  }
  return(zscore_matrix)
}

# Calculate the Z-scores based on enrichment factors
go_zscores <- calculate_go_zscores_from_enrichment(go_results)

# Remove rows (GO terms) with all Z-scores as NA
go_zscores <- go_zscores[apply(go_zscores, 1, function(row) any(!is.na(row))), ]

# Optional: Filter to keep top N GO terms based on variance across clusters
n_top_terms <- min(50, nrow(go_zscores)) # Adjust the number as needed
if (nrow(go_zscores) > 1) {
  term_variance <- apply(go_zscores, 1, var, na.rm = TRUE)
  top_n_terms <- names(sort(term_variance, decreasing = TRUE)[1:n_top_terms])
  go_zscores_to_plot <- go_zscores[top_n_terms, ]
} else {
  go_zscores_to_plot <- go_zscores
}



# STEP 8: Combine GO enrichment results into one data frame
go_combined <- bind_rows(go_results)

# STEP 9: Create a matrix of enrichment scores (-log10 adjusted p-values)
go_matrix <- go_combined %>%
  mutate(Zscore = -log10(p.adjust)) %>%
  select(cluster, Description, Zscore) %>%
  pivot_wider(names_from = cluster, values_from = Zscore, values_fill = list(Zscore = 0))

# Select top 30 GO terms by average Z-score across clusters
go_matrix <- go_matrix %>%
  mutate(mean_z = rowMeans(select(., -Description))) %>%
  arrange(desc(mean_z)) %>%
  slice_head(n = 30) %>%
  select(-mean_z)

# STEP 10: Convert the dataframe to a matrix for heatmap plotting
go_mat <- as.matrix(go_matrix %>% select(-Description))
rownames(go_mat) <- go_matrix$Description

# STEP 11: Define a color scale for the heatmap
col_fun <- colorRamp2(
  c(0, 25, 50),
  c("white", "red", "darkred")
)

# STEP 12: Plot the heatmap and save it to a PDF file
pdf("Plots/subclustered/GOHeatmap_TopMarkers.pdf", width = 12)
Heatmap(
  go_mat,
  name = "Enrichment score\n(-log10 adjusted p-value)",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title_position = "topcenter",
    legend_direction = "horizontal"
  ),
  border = TRUE,
  column_title = "GO Enrichment of Top Markers per Cluster",
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()





# Prepare the data by counting cells per timepoint and batch
cell_counts <- full_dataset@meta.data %>%
  group_by(timepoint, batch) %>%
  summarise(cell_count = n()) %>%
  ungroup()

# Plot 5: Bar plot showing number of cells per timepoint, split by batch
plot5 <- ggplot(cell_counts, aes(x = timepoint, y = cell_count, fill = batch)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Number of Cells per Timepoint",
       subtitle = "Stacked by Batch",
       x = "Timepoint",
       y = "Number of Cells",
       fill = "Batch") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"), 
        axis.title.y = element_text(size = 14, face = "bold"),  
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) 
# Save the bar plot
ggsave(filename = paste0(output_folder, "/5_cells_per_timepoint_split_by_batch.pdf"), plot = plot5, width = 10, height = 6)

# Define the marker genes for each group
marker_genes <- list(
  "Neuronal Progenitors" = c("dpn", "tll"),
  "Sensory Progenitors" = c("pros", "ase"),
  "Glia" = c("repo", "gcm"),
  "Cholinergic Neurons" = c("VAChT", "ChAT"),
  "Glutamatergic Neurons" = c("VGlut"),
  # "Peptidergic Neurons" = c("FMRFa", "ITP", "AstC-R2"),
  "Monoaminergic Neurons" = c("Tdc2", "ple", "Vmat", "Trh", "Hdc", "Tbh")
)

# Combine the markers into a single vector
all_marker_genes <- unlist(marker_genes)

# Set the Seurat object identity to "seurat_clusters"
Idents(full_dataset) <- full_dataset@meta.data$seurat_clusters

# Generate the Dot Plot in Seurat with rotated x-axis labels
dot_plot <- DotPlot(full_dataset, features = all_marker_genes, scale.max = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the plot with adjusted size and resolution
ggsave(filename = paste0(output_folder, "/6_dotplot_marker_genes.pdf"), plot = dot_plot, width = 10, height = 20)

# Set the Seurat object identity to "seurat_clusters"
Idents(full_dataset) <- full_dataset@meta.data$bdgp_putatative_celltype

# Generate the Dot Plot in Seurat with rotated x-axis labels
dot_plot <- DotPlot(full_dataset, features = all_marker_genes, scale.max = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the plot with adjusted size and resolution
ggsave(filename = paste0(output_folder, "/6A_dotplot_marker_genes.pdf"), plot = dot_plot, width = 10, height = 20)

monoaminergic_markers <- list(
  "Serotonergic" = c("Vmat", "Trh", "Ddc", "SerT"),
  "Dopaminergic" = c("ple", "DAT"),
  "Histaminergic" = c("Hdc"),
  "Tyraminergic" = c("Tdc1", "Tdc2"),
  "Octopaminergic" = c("Tbh")
)

monoaminergic_receptors <- list(
  "Serotonergic" = c("5-HT1A", "5-HT1B", "5-HT2A", "5-HT2B", "5-HT7"),
  "Dopaminergic" = c("Dop1R1", "Dop1R2", "Dop2R", "DopEcR"),
  "Histaminergic" = c("HisCl1"),
  "Tyraminergic" = c("TyrR", "TyrRII"),
  "Octopaminergic" = c("Oct-TyrR", "Octalpha2R", "Octbeta1R", "Octbeta2R", "Octbeta3R")
)

# Combine the markers into a single vector
all_monoaminergic_marker_genes <- unlist(monoaminergic_markers)

# Combine the markers into a single vector
all_monoaminergic_markers <- unlist(c(monoaminergic_markers, monoaminergic_receptors))
all_monoaminergic_marker_genes <- c("Vmat", "Trh", "Ddc", "SerT", "ple", "DAT", "Hdc", "Tdc1", "Tdc2", "Tbh", "5-HT1A", "5-HT1B", "5-HT2A", "5-HT2B", "5-HT7", "Dop1R1", "Dop1R2", "Dop2R", "DopEcR", "HisCl1", "TyrR", "TyrRII", "Oct-TyrR", "Octalpha2R", "Octbeta1R", "Octbeta2R", "Octbeta3R")

# Set the Seurat object identity to "seurat_clusters"
Idents(full_dataset) <- full_dataset@meta.data$seurat_clusters

# Generate the Dot Plot in Seurat with rotated x-axis labels
dot_plot <- DotPlot(full_dataset, features = all_monoaminergic_marker_genes, scale.max = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the plot with adjusted size and resolution
ggsave(filename = paste0(output_folder, "/7_dotplot_monoaminergic_marker_genes.pdf"), plot = dot_plot, width = 10, height = 20)

qsave(full_dataset, "/DataDrives/Drive3/Clifton/FullDrosophilaProject/data/full_analysis/full_dataset.qs")
use_condaenv("sceasy")
sceasy::convertFormat(obj = full_dataset, from = "seurat", to = "anndata", outFile = "/DataDrives/Drive3/Clifton/FullDrosophilaProject/data/full_analysis/full_dataset.h5ad")

########  GO Enrichment Analysis ######## 
# Load required libraries
library(Seurat)
library(clusterProfiler)
library(org.Dm.eg.db) 
library(dplyr)
library(ggplot2)
library(enrichplot)
library(cowplot)  # For arranging multiple plots in a grid layout

# Step 1: Identify marker genes for each cluster
Idents(full_dataset) <- full_dataset@meta.data$bdgp_putatative_celltype_unique
markers <- FindAllMarkers(full_dataset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top markers
head(markers)
qsave(markers, "data/bdgp_putatative_celltype_unique_findallmarkers.qs")

# Convert gene symbols to Entrez IDs using the org.Dm.eg.db database for Drosophila
# Select top markers for each cluster
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(top_markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)

# Merge the Entrez IDs with the original markers data
top_markers_entrez <- merge(top_markers, gene_entrez, by.x = "gene", by.y = "SYMBOL")

# Step 2: Perform GO enrichment analysis for each cluster using clusterProfiler
go_results_list <- list()

for (cluster in unique(top_markers_entrez$cluster)) {
  
  # Subset markers for this specific cluster
  cluster_genes <- top_markers_entrez %>% filter(cluster == !!cluster) %>% pull(ENTREZID)
  
  # Check if there are sufficient genes to perform GO analysis
  if (length(cluster_genes) > 0) {
    
    # Run GO enrichment analysis
    go_enrichment <- enrichGO(gene = cluster_genes, 
                              OrgDb = org.Dm.eg.db, 
                              keyType = "ENTREZID", 
                              ont = "BP",   # Biological process (BP), can also be MF or CC
                              pAdjustMethod = "BH", 
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.05)
    
    # Only store results if there are significant GO terms found
    if (!is.null(go_enrichment) && nrow(go_enrichment) > 0) {
      go_results_list[[cluster]] <- go_enrichment
      
      # Sanitize the cluster name to avoid issues with special characters in file names
      sanitized_cluster <- gsub("[^[:alnum:]_]", "_", cluster)
      
      # Save results to a CSV file for each cluster
      write.csv(as.data.frame(go_enrichment), file = paste0("GO_results_cluster_", sanitized_cluster, ".csv"))
      
    } else {
      message(paste("No GO terms found for cluster:", cluster))
    }
  }
}

# Step 3: Create and plot GO enrichment bar plots for each cluster
# Function to generate a bar plot for the top GO terms of a cluster
create_go_barplot <- function(go_result, cluster) {
  if (!is.null(go_result) && nrow(go_result) > 0) {
    p <- barplot(go_result, showCategory = 10, title = paste("Cluster", cluster)) +
      theme_classic() +
      theme(plot.title = element_text(size = 16, face = "bold"),
            axis.text.x = element_text(size = 14))
    return(p)
  } else {
    return(NULL)
  }
}

# Step 4: Create a list of bar plots for each cluster
plot_list <- list()

for (cluster in unique(top_markers_entrez$cluster)) {
  if (!is.null(go_results_list[[cluster]])) {
    go_barplot <- create_go_barplot(go_results_list[[cluster]], cluster)
    if (!is.null(go_barplot)) {
      plot_list[[cluster]] <- go_barplot
    }
  }
}

# Step 5: Combine the individual bar plots into a grid layout
# Adjust ncol based on the number of clusters
final_plot <- plot_grid(plotlist = plot_list, ncol = 8)  # You can change ncol for different layouts

# Step 6: Save the final plot
ggsave("plots/full_dataset_analysis/GO_enrichment_clusters_combined.pdf", plot = final_plot, width = 48, height = 30)

# Convert gene symbols to Entrez IDs using the org.Dm.eg.db database for Drosophila
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
gene_entrez <- bitr(top_markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)
top_markers_entrez <- merge(top_markers, gene_entrez, by.x = "gene", by.y = "SYMBOL")

# Step 2: Prepare a list of genes for each cluster
gene_list <- split(top_markers_entrez$ENTREZID, top_markers_entrez$cluster)

# Step 3: Perform GO enrichment analysis for all clusters together using compareCluster
go_enrichment_clusters <- compareCluster(
  geneCluster = gene_list, 
  fun = "enrichGO", 
  OrgDb = org.Dm.eg.db, 
  ont = "BP",  # Biological process
  keyType = "ENTREZID", 
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Step 4: Create a combined dot plot for GO enrichment across all clusters
go_dotplot <- dotplot(go_enrichment_clusters, showCategory = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 5: Save the GO enrichment plot
ggsave("plots/full_dataset_analysis/GO_enrichment_all_clusters_dotplot.pdf", plot = go_dotplot, width = 15, height = 10)