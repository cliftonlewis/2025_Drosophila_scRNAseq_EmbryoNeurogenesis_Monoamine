# ============================
# Integrated Analysis Script
# ============================

# ────────────────────────────────────────────────────
# 1. Environment setup
# ────────────────────────────────────────────────────
.libPaths('~/R/x86_64-pc-linux-gnu-library/4.5/')
options(future.globals.maxSize = 40000 * 1024^2)
setwd("/DataDrives/Drive2/Clifton/R_Projects/2025_Drosophila_scRNAseq_MonoamineSpecification/ANALYSIS/Step_2_Initial_preprocessing/Step_2.4_Mature_neuronal_subset_preprocessing/")

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

# Load datasets
mature_neuronal_subset <- qread(file = "Data/Timepoint_Datasets/mature_neuronal_subset/1_mature_neuronal_subset_preannotated.qs")
sample_name <- "mature_neuronal_subset"
pcs <- 129

mature_neuronal_subset <- run_normalisation(mature_neuronal_subset)
mature_neuronal_subset <- cell_cycle_scoring(sample_object = mature_neuronal_subset)
mature_neuronal_subset <- generate_variable_features_plot(mature_neuronal_subset, sample_name)

# Run PCA
mature_neuronal_subset <- RunPCA(object = mature_neuronal_subset, features = VariableFeatures(object = mature_neuronal_subset), npcs = 200)

# Plot the elbow plot
plot1 <- ElbowPlot(object = mature_neuronal_subset, ndims=200)
print(plot1)
ggsave(paste0("Plots/Timepoint_Datasets/", sample_name, "/4_", sample_name,  "_ElbowPlot.png"), plot1, width = 14, height = 6)

# Calculate the percentage of variation associated with each PC
pct <- mature_neuronal_subset[["pca"]]@stdev / sum(mature_neuronal_subset[["pca"]]@stdev) * 100

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

# Run the JackStraw procedure for determining statistically significant PCs
mature_neuronal_subset <- JackStraw(mature_neuronal_subset, dims = 200)
mature_neuronal_subset <- ScoreJackStraw(mature_neuronal_subset, dims = 1:200)

# JackStraw plot to visualize the statistical significance of each PC
pcs_jackstraw <- length(which(mature_neuronal_subset[["pca"]]@jackstraw@overall.p.values[, "Score"] < 0.05))

plot2 <- JackStrawPlot(mature_neuronal_subset, dims = 1:200) +
  annotate("text", x = pcs_jackstraw + 3, y = max(mature_neuronal_subset[["pca"]]@jackstraw@overall.p.values[, "Score"]), 
           label = paste("PCs chosen: ", pcs_jackstraw), color = "blue", size = 5) +  # Annotate chosen PCs
  ggtitle("JackStraw Plot: Significance of Principal Components") +
  labs(x = "Principal Components", y = "p-value")
print(plot2)
# Save the JackStraw plot
ggsave(paste0("Plots/Timepoint_Datasets/", sample_name, "/5B_", sample_name, "_JackStrawPlot_significant.png"), plot2, width = 20, height = 6)

print(paste("Optimal number of PCs based on JackStraw plot:", pcs_jackstraw))

# Find neighbors
mature_neuronal_subset <- FindNeighbors(mature_neuronal_subset, dims = 1:pcs)

# Perform clustering
resolution.range <- seq(from = 0, to = 5, by = 0.1)
mature_neuronal_subset <- FindClusters(mature_neuronal_subset,
                                       reduction = "pca", dims = 1:pcs,
                                       resolution = resolution.range,
                                       save.SNN = TRUE,
                                       print.output = FALSE)
# Extract metadata columns corresponding to clustering resolutions
cluster_cols <- grep("^RNA_snn_res\\.", colnames(mature_neuronal_subset@meta.data), value = TRUE)

# Summarize number of clusters per resolution
cluster_summary <- sapply(cluster_cols, function(x) length(unique(mature_neuronal_subset@meta.data[[x]])))

# Put into a data frame
res_cluster_df <- data.frame(
  resolution = gsub("RNA_snn_res.", "", cluster_cols),
  n_clusters = cluster_summary
)

print(res_cluster_df)

# Generate clustree plot
plot1 <- clustree(mature_neuronal_subset) + theme_light()
ggsave(paste0("Plots/Timepoint_Datasets/", sample_name, "/6_", sample_name, "_Clustree.png"), plot1, width = 15, height = 15)

# Generate clustree plot with SC3 stability color-coding
plot2 <- clustree(mature_neuronal_subset, node_colour = "sc3_stability")
ggsave(paste0("Plots/Timepoint_Datasets/", sample_name, "/7_", sample_name, "_Clustree_sc3_stability.png"), plot2, width = 15, height = 15)

mature_neuronal_subset <- FindNeighbors(mature_neuronal_subset, dims = 1:pcs)
mature_neuronal_subset <- FindClusters(mature_neuronal_subset, reduction = "pca", dims = 1:pcs, resolution = 1.9)
mature_neuronal_subset <- run_dimred(mature_neuronal_subset, sample_name, pcs)
qsave(mature_neuronal_subset, file = paste0("Data/Timepoint_Datasets/", sample_name, "/2_", sample_name, "_preprocessed.qs"))

# Verify the merged object
print(mature_neuronal_subset)
head(mature_neuronal_subset@meta.data)

mature_neuronal_subset <- qread(paste0("Data/Timepoint_Datasets/", sample_name, "/2_", sample_name, "_preprocessed.qs"))
mature_neuronal_subset <- mature_neuronal_subset 

metadata <- mature_neuronal_subset@meta.data

# Function to get the most frequent neuronal_annotation_broad for each seurat_cluster
get_most_frequent_label <- function(cluster_id) {
  # Subset the metadata for the specific cluster
  cluster_data <- metadata[metadata$seurat_clusters == cluster_id, ]
  
  # Filter out NA or NULL values in neuronal_annotation_broad
  cluster_data <- cluster_data[!is.na(cluster_data$neuronal_annotation_broad) & cluster_data$neuronal_annotation_broad != "NULL", ]
  
  # Check if we have any valid values left after filtering
  if (nrow(cluster_data) == 0) {
    return("Unknown")  # Assign "Unknown" if no valid cell type exists for the cluster
  }
  
  # Find the most frequent neuronal_annotation_broad within this cluster
  most_frequent <- names(sort(table(cluster_data$neuronal_annotation_broad), decreasing = TRUE))[1]
  
  return(most_frequent)
}

# Apply the function to each unique seurat_cluster
most_frequent_labels <- sapply(unique(metadata$seurat_clusters), get_most_frequent_label)

# Create a mapping of seurat_cluster -> most frequent neuronal_annotation_broad
cluster_to_celltype <- setNames(most_frequent_labels, unique(metadata$seurat_clusters))

# Now create a new column combining seurat_cluster and most frequent neuronal_annotation_broad
metadata$neuronal_annotation_broad <- paste0(cluster_to_celltype[as.character(metadata$seurat_clusters)])

# Now create a new column combining seurat_cluster and most frequent celltype
metadata$neuronal_annotation_broad_unique <- paste0(
  metadata$seurat_clusters, "_", cluster_to_celltype[as.character(metadata$seurat_clusters)]
)

# **Corrected Assignment:** Add the updated metadata back to the Seurat object
mature_neuronal_subset@meta.data$neuronal_annotation_broad <- metadata$neuronal_annotation_broad
mature_neuronal_subset@meta.data$neuronal_annotation_broad_unique <- metadata$neuronal_annotation_broad_unique

# Optionally, view the result
table(mature_neuronal_subset@meta.data$neuronal_annotation_broad)
table(mature_neuronal_subset@meta.data$neuronal_annotation_broad_unique)

# Plot 3: Annotated by Timepoint
Idents(mature_neuronal_subset) <- mature_neuronal_subset$timepoint

# Load your cluster-to-color mapping from a CSV file
# The CSV file should have two columns: "Label" and "Color"
annotation_table <- read.csv("../Supplementary_Data/timepoint2colour_mapping.csv")

# Ensure the Seurat object has cluster identities assigned
clusters_in_seurat <- levels(Idents(mature_neuronal_subset))  # Get unique clusters in the Seurat object

# Subset the annotation table to include only clusters present in the Seurat object
filtered_annotation_table <- annotation_table[annotation_table$Label %in% clusters_in_seurat, ]

# Create a named vector for colors based on the filtered annotation table
color_mapping_timepoint <- setNames(filtered_annotation_table$Color, filtered_annotation_table$Label)

# Map colors to each cell in the Seurat object
cell_colors <- color_mapping_timepoint[as.character(mature_neuronal_subset@meta.data$timepoint)]

# Add the colors as metadata in the Seurat object
mature_neuronal_subset@meta.data$timepoint_color <- cell_colors

# Verify the new metadata column
head(mature_neuronal_subset@meta.data)

plot3 <- dittoDimPlot(object = mature_neuronal_subset, 'timepoint', reduction.use = "umap", color.panel = color_mapping_timepoint) +
  theme_minimal() +
  labs(title = "UMAP Plot", subtitle = "Neuronal Dataset Annotated by Timepoint") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "right")
print(plot3)
ggsave(filename = paste0("Plots/Timepoint_Datasets/mature_neuronal_subset/",   "/3_mature_neuronal_dataset_timepoint_UMAP2D.pdf"), plot = plot3, width = 10, height = 6)

# Plot 4: Annotated 
Idents(mature_neuronal_subset) <- mature_neuronal_subset$neuronal_annotation_broad
plot4 <- dittoDimPlot(object = mature_neuronal_subset, 'neuronal_annotation_broad', reduction.use = "umap", do.label = TRUE) +
  theme_minimal() +
  labs(title = "UMAP Plot", subtitle = "Neuronal Dataset Annotated") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "right")
print(plot4)
ggsave(filename = paste0("Plots/Timepoint_Datasets/mature_neuronal_subset/",   "/4_mature_neuronal_dataset_neuronal_annotation_broad_UMAP2D.pdf"), plot = plot4, width = 10, height = 6)


Idents(mature_neuronal_subset) <- mature_neuronal_subset$seurat_clusters
plot4 <- dittoDimPlot(object = mature_neuronal_subset, 'seurat_clusters', reduction.use = "umap", do.label = TRUE) +
  theme_minimal() +
  labs(title = "UMAP Plot", subtitle = "Neuronal Dataset - Clusters") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "right")
print(plot4)
ggsave(filename = paste0("Plots/Timepoint_Datasets/mature_neuronal_subset/",   "/4_mature_neuronal_dataset_neuronal_annotation_broad_UMAP2D.pdf"), plot = plot4, width = 10, height = 6)


qsave(mature_neuronal_subset, file = "Data/Timepoint_Datasets/mature_neuronal_subset/3_mature_neuronal_subset_annotated_w_unique.qs")

## Neuronal lineage
mature_neuronal_subset <- qread("Data/Timepoint_Datasets/mature_neuronal_subset/3_mature_neuronal_subset_annotated_w_unique.qs")

# Define your grouped features as a named list
nt_df <- list(
  "Neuroblasts" = c("trbl", "dpn", "wor", "ase", "mira", "insc", "stg", "tll"),
  "GMCs" = c("tap", "dap", "pros", "numb"),
  "New-born neurons" = c("Hey", "E(spl)m6-BFM"),
  "Immature neurons" = c("elav", "CadN", "fne"),
  "Mature neurons" = c("brp", "nSyb", "Syn", "Syt1", "Syt4"),
  "Acetylcholine" = c("ChAT", "VAChT", "ChT"),
  "Glutamate" = c("VGlut"),
  "GABA" = c("Gad1", "VGAT", "Gat"),
  "Serotonin" = c("Trh", "SerT"),
  "Dopamine" = c("ple", "DAT"),
  "Monoamine" = c("Tbh", "Tdc1", "Tdc2", "Ddc", "Vmat")
  #, "AcetylcholineR" = c(
  #   "nAChRalpha1",
  #   "nAChRalpha2",
  #   "nAChRalpha3",
  #   "nAChRalpha4",
  #   "nAChRalpha5",
  #   "nAChRalpha6",
  #   "nAChRalpha7",
  #   "nAChRbeta1",
  #   "nAChRbeta2",
  #   "nAChRbeta3",
  #   "mAChR-A",
  #   "mAChR-B",
  #   "mAChR-C"
  # ),
  # "GlutamateR" = c("GluRIIA", "GluRIIB", "GluRIIC", "GluRIID", "GluRIIE", "mGluR"),
  # "GABAR" = c("Rdl", "Lcch3", "GABA-B-R1", "GABA-B-R2", "GABA-B-R3"),
  # "SerotoninR" = c("5-HT1A", "5-HT1B", "5-HT2A", "5-HT2B", "5-HT7"),
  # "DopamineR" = c("Dop1R1", "Dop1R2", "Dop2R", "DopEcR"),
  # "HistamineR" = c("HisCl1"),
  # "TyramineR" = c("TyrR", "TyrRII"),
  # "OctopamineR" = c(
  #   "Oct-TyrR",
  #   "Octalpha2R",
  #   "Octbeta1R",
  #   "Octbeta2R",
  #   "Octbeta3R"
  # )
)

plot5 <- DotPlot2(seu = mature_neuronal_subset, features = nt_df, group.by = "seurat_clusters", color_scheme = "BuRd")
print(plot5)
ggsave(filename = paste0("Plots/Timepoint_Datasets/mature_neuronal_subset/",   "/5_mature_neuronal_dataset_neuro_markers_seurat_clusters.pdf"), plot = plot5, width = 20, height = 20)

plot5 <- DotPlot2(seu = mature_neuronal_subset, features = nt_df, group.by = "celltype_merge", color_scheme = "BuRd")
print(plot5)
ggsave(filename = paste0("Plots/Timepoint_Datasets/mature_neuronal_subset/",   "/5_mature_neuronal_dataset_neuro_markers_celltype_merge.pdf"), plot = plot5, width = 6, height = 25)

plot5 <- DotPlot2(seu = mature_neuronal_subset, features = nt_df, group.by = "neuronal_annotation_broad", color_scheme = "BuRd")
print(plot5)
ggsave(filename = paste0("Plots/Timepoint_Datasets/mature_neuronal_subset/",   "/5_mature_neuronal_dataset_neuro_markers_neuronal_annotation_broad.pdf"), plot = plot5, width = 6, height = 20)

plot5 <- DotPlot2(seu = mature_neuronal_subset, features = nt_df, group.by = "timepoint", color_scheme = "BuRd")
print(plot5)
ggsave(filename = paste0("Plots/Timepoint_Datasets/mature_neuronal_subset/",   "/5_mature_neuronal_dataset_timepoint.pdf"), plot = plot5, width = 10, height = 20)

# New neuronal annotations

# Ensure the active identities are set to the clusters
Idents(mature_neuronal_subset) <- mature_neuronal_subset$seurat_clusters

# Check the current levels to verify the cluster names
levels(mature_neuronal_subset)

# Define the mapping vector for new cluster annotations
new.cluster.ids <- c(
  "0" = "Neuroblasts/GMCs/Immature_neurons",
  "1" = "Neuroblasts/GMCs",
  "2" = "Acetylcholine",
  "3" = "GMCs/New-born_neurons/Immature_neurons",
  "4" = "Acetylcholine/GABA",
  "5" = "Neuroblasts/GMCs/Immature_neurons",
  "6" = "Acetylcholine/GABA",
  "7" = "Acetylcholine/GABA",
  "8" = "GABA",
  "9" = "Neuroblasts/GMCs",
  "10" = "Neuroblasts/GMCs/Immature_neurons",
  "11" = "Neuroblasts/GMCs/Immature_neurons",
  "12" = "Immature_neurons",
  "13" = "Neuroblasts",
  "14" = "GABA/Glutamate",
  "15" = "Neuroblasts/GMCs/Immature_neurons",
  "16" = "Neuroblasts/GMCs/Immature_neurons",
  "17" = "Neuroblasts/GMCs",
  "18" = "Neuroblasts/GMCs",
  "19" = "Neuroblasts/GMCs",
  "20" = "Neuroblasts/GMCs",
  "21" = "Neuroblasts",
  "22" = "GMCs/New-born_neurons/Immature_neurons",
  "23" = "Neuroblasts/GMCs/New-born_neurons/Immature_neurons",
  "24" = "Acetylcholine/GABA",
  "25" = "Neuroblasts/GMCs",
  "26" = "Neuroblasts",
  "27" = "Neuroblasts/GMCs",
  "28" = "Unknown",
  "29" = "New-born_neurons/Immature_neurons",
  "30" = "Glutamate",
  "31" = "Unknown",
  "32" = "Neuroblasts",
  "33" = "Neuroblasts", 
  "34" = "Neuroblasts/GMCs/New-born_neurons/Immature_neurons",
  "35" = "Monoamine",  
  "36" = "Monoamine/GABA",
  "37" = "Neuroblasts/GMCs/New-born_neurons/Immature_neurons",
  "38" = "Acetylcholine",
  "39" = "Neuroblasts/GMCs/Immature_neurons",
  "40" = "Neuroblasts/GMCs/Immature_neurons",
  "41" = "Neuroblasts/GMCs",
  "42" = "Neuroblasts/GMCs/Immature_neurons",
  "43" = "Acetylcholine"
)

mature_neuronal_subset@meta.data$neuronal_annotation_fine <- mature_neuronal_subset@meta.data$seurat_clusters
Idents(mature_neuronal_subset) <- mature_neuronal_subset$neuronal_annotation_fine
levels(mature_neuronal_subset)

# Rename the active identities in the Seurat object
mature_neuronal_subset <- RenameIdents(mature_neuronal_subset, new.cluster.ids)

# Save the new identities into a metadata column for future reference.
mature_neuronal_subset@meta.data$neuronal_annotation_fine <- Idents(mature_neuronal_subset)
table(mature_neuronal_subset@meta.data$neuronal_annotation_fine)
# Verify the changes by checking the metadata
head(mature_neuronal_subset@meta.data)
shared_cell_types <- c(
  "Acetylcholine",
  "Acetylcholine/GABA",
  "GABA",
  "GABA/Glutamate",
  "GABA/Serotonin",
  "Glutamate",
  "Serotonin",
  "Serotonin/GABA",
  "Monoamine",
  "Monoamine/Acetylcholine",
  "Monoamine/Serotonin",
  "Monoamine/GABA",
  "Immature_neurons",
  "New-born_neurons/Immature_neurons",
  "Unknown_mature_neurons",
  "Neuroblasts",
  "Neuroblasts/GMCs",
  "Neuroblasts/GMCs/Immature_neurons",
  "GMCs",
  "GMCs/New-born_neurons/Immature_neurons",
  "Neuroblasts/GMCs/New-born_neurons/Immature_neurons",
  "Unknown"
)

shared_color_palette <- c(
  # Acetylcholine-related
  "Acetylcholine"                        = "#FFD700",  # vivid gold
  "Acetylcholine/GABA"                   = "#FFEE58",  # sunflower yellow
  
  # GABA-related (varied reds and pinks)
  "GABA"                                 = "#B71C1C",  # dark red
  "GABA/Glutamate"                       = "#D84315",  # burnt orange-red
  "GABA/Serotonin"                       = "#F06292",  # deep pink
  "Serotonin/GABA"                       = "#E91E63",  # strong pink-rose
  "Monoamine/GABA"                       = "#EF5350",  # soft red-pink
  
  # Glutamate-related
  "Glutamate"                            = "#43A047",  # strong green
  
  # Serotonin (purple family)
  "Serotonin"                            = "#8E24AA",  # deep purple
  
  # Monoamine-related (orange shades)
  "Monoamine"                            = "#FB8C00",  # bright orange
  "Monoamine/Acetylcholine"              = "#FFA726",  # soft orange
  "Monoamine/Serotonin"                  = "#FF7043",  # orange-coral
  
  # Neuroblasts & GMCs (cool blues)
  "Neuroblasts"                          = "#1565C0",  # cobalt blue
  "Neuroblasts/GMCs"                     = "#1E88E5",  # vivid blue
  "Neuroblasts/GMCs/Immature_neurons"    = "#64B5F6",  # sky blue
  "GMCs"                                 = "#0D47A1",  # navy
  "GMCs/New-born_neurons/Immature_neurons" = "#1976D2",  # medium blue
  "Neuroblasts/GMCs/New-born_neurons/Immature_neurons" = "#90CAF9",  # pale blue
  
  # Developmental/Immature/Unknown
  "Immature_neurons"                     = "#29B6F6",  # bright cyan
  "New-born_neurons/Immature_neurons"    = "#4DD0E1",  # teal
  "Unknown_mature_neurons"               = "#757575",  # neutral gray
  "Unknown"                              = "#BDBDBD"   # light gray
)

plot7 <- DimPlot2(
  mature_neuronal_subset,
  features = c("neuronal_annotation_fine"),
  cols = list("neuronal_annotation_fine" = shared_color_palette), pt.size = 1.25) +
  labs(title = "UMAP Plot") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    legend.position = "right"
  )

print(plot7)
ggsave(filename = paste0("Plots/Timepoint_Datasets/mature_neuronal_subset/",   "/7_neuronal_dataset_neuronal_annotation_fine_UMAP2D.pdf"), plot = plot7, width = 10, height = 6)

qsave(mature_neuronal_subset, file = "Data/Timepoint_Datasets/mature_neuronal_subset/4_mature_neuronal_subset_annotated.qs")


mature_neuronal_subset <- qread("Data/Timepoint_Datasets/mature_neuronal_subset/4_mature_neuronal_subset_annotated.qs")
table(mature_neuronal_subset@meta.data$neuronal_annotation_fine)
library(reticulate)
use_condaenv('sccustomize')
as.anndata(mature_neuronal_subset, file_path = "Data/Timepoint_Datasets/mature_neuronal_subset/", "4_mature_neuronal_subset_annotated.h5ad")
