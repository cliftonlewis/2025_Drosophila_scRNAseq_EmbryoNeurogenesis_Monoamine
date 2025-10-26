###########################################
# INSTALLATION & ENVIRONMENT SETUP
###########################################
# Set the library path for custom R packages and increase the size limit for future globals.
.libPaths('~/R/x86_64-pc-linux-gnu-library/4.5/')
options(future.globals.maxSize = 40000 * 1024^2)

# Set the working directory.
setwd("/DataDrives/Drive2/Clifton/R_Projects/2025_Drosophila_scRNAseq_MonoamineSpecification/ANALYSIS/Step_2_Initial_preprocessing/Step_2.1_Sample_preprocessing/")

# Define directories for saving processed data and generated plots.
data_dir <- "Data/Timepoint_Datasets"
plots_dir <- "Plots/Timepoint_Datasets"

# Create base directories if they don't already exist.
dir.create("Data", showWarnings = FALSE)
dir.create("Plots", showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

###########################################
# LOAD REQUIRED PACKAGES
###########################################
# Load scRNA-seq analysis packages.
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
# Uncomment if using sceasy for conversions.
# library(sceasy)
library(org.Dm.eg.db)
library(biomaRt)
library(RAPToR)
library(drosoRef)

# Load data science and utility packages.
library(tidyverse)
library(readr)
library(stringr)
library(dplyr)
library(qs)
library(magrittr)

# Load packages for plotting and aesthetics.
library(patchwork)
library(ggplot2)
library(cowplot)
library(dittoSeq)

# Load additional packages for annotation and citation management.
organism <- "org.Dm.eg.db"
library(organism, character.only = TRUE)
library(tibble)
library(tools)  # For converting citation information to BibTeX format.

###########################################
# EXPORT PACKAGE CITATIONS
###########################################
# This function extracts citation information for all loaded packages.
get_package_citations <- function() {
  pkg_info <- sessionInfo()$otherPkgs
  citations_df <- tibble::tibble(
    Package = character(),
    Version = character(),
    Citation = character()
  )
  # Loop through each package and try to retrieve its citation.
  for (pkg in names(pkg_info)) {
    version <- pkg_info[[pkg]]$Version
    citation_text <- tryCatch({
      citation_info <- citation(pkg)
      bibtex <- toBibtex(citation_info)
      paste(bibtex, collapse = "\n")
    }, error = function(e) NA_character_)
    citations_df <- rbind(citations_df, tibble::tibble(
      Package = pkg,
      Version = version,
      Citation = citation_text
    ))
  }
  return(citations_df)
}

# Save package citations to a CSV file for reference.
citation_data <- get_package_citations()
print(citation_data)
write.csv(citation_data, '../Citations/citations.csv', row.names = FALSE, quote = FALSE)

###########################################
# DEFINE SAMPLE METADATA & SETUP DIRECTORIES
###########################################
# Define vectors with sample IDs, timepoints, and batch information.
sample <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8",
            "A9", "A10", "A11", "A12", "A13", "A14", "A15", "A16")
timepoint <- c("hrs_00_03", "hrs_02_05", "hrs_04_07", "hrs_05_08",
               "hrs_06_10", "hrs_09_13", "hrs_12_17", "hrs_16_22",
               "hrs_00_03", "hrs_02_05", "hrs_04_07", "hrs_05_08",
               "hrs_06_10", "hrs_09_13", "hrs_12_17", "hrs_16_22")
# Use replicates for batches (each replicate covers 8 samples).
batch <- rep(c("Replicate 1", "Replicate 2"), each = 8)

# Create a metadata dataframe to store sample-level information.
metadata_df <- data.frame(
  sample = sample,
  timepoint = timepoint,
  batch = batch,
  stringsAsFactors = FALSE
)

# Specify the sample name for this analysis (example: sample A1).
sample_name <- "hrs_00_03_A9_Rep2"
# Create directories specific to the sample for data and plots.
sample_dir_data <- file.path(data_dir, sample_name)
sample_dir_plots <- file.path(plots_dir, sample_name)
dir.create(sample_dir_data, showWarnings = FALSE, recursive = TRUE)
dir.create(sample_dir_plots, showWarnings = FALSE, recursive = TRUE)
file_path2 <- paste0("/DataDrives/Drive2/Clifton/R_Projects/2025_Drosophila_scRNAseq_MonoamineSpecification/ANALYSIS/Step_2_Initial_preprocessing")

###########################################
# SOURCE CUSTOM FUNCTIONS
###########################################
# Load custom functions required for data processing.

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
source(file.path(file_path2, "Functions/doubletfinder.R"))
source(file.path(file_path2, "Functions/fishers_test.R"))  

###########################################
# SAMPLE ANALYSIS STEPS
###########################################
# Set additional parameters and file paths.
sample_id <- "A9"
chosen_resolution <- 2.7
file_path <- paste0("/DataDrives/Drive2/Clifton/R_Projects/2025_Drosophila_scRNAseq_MonoamineSpecification/ANALYSIS/Step_1_CellRanger_and_Velocyto_Processing/", sample_name)
bdgp_db <- read.csv(file.path(file_path2, "Annotation_Tables/insitu_annot_hrs_00_03.csv"), header = TRUE)
is_loom <- TRUE  # Indicate whether the data is in loom format or not.
low_counts_threshold <- 800
high_counts_threshold <- 100000
low_genes_threshold <- 200
high_genes_threshold <- 8000
mt_threshold <- 8  # Mitochondrial gene threshold (%)
exp_multiplet <- 0.039  # Expected multiplet rate
pcs <- 38

# Ensure metadata is provided before proceeding.
if (is.null(metadata_df)) stop("Metadata dataframe 'metadata_df' is required.")

# Set random seed for reproducibility.
set.seed(123)

# --- Step 1: Load Data into Seurat Object ---
cat("Step 1: Loading data into Seurat object\n")
if (is_loom) {
  sample_object <- create_seurat_object_loom(sample_name, file_path)
} else {
  sample_object <- create_seurat_object_standard(sample_name, sample_id, file_path, file_path2 = file_path2)
}
cat("Seurat object created successfully.\n")

# --- Step 2: Add Metadata ---
cat("Step 2: Adding metadata to Seurat object\n")
sample_object <- add_metadata(sample_object, metadata_df, sample_id, file_path2 = file_path2)
cat("Metadata added successfully.\n")

# --- Step 3: Perform Quality Control (QC) and Generate QC Plots ---
cat("Step 3: Performing QC and generating plots\n")
sample_object <- perform_qc_and_plots(sample_object, sample_name, file_path2 = file_path2)
cat("QC and plots generated successfully.\n")

# --- Step 4: Filter Cells Based on QC Metrics ---
cat("Step 4: Performing QC filtering\n")
sample_object <- perform_qc_filtering(sample_object, sample_name, 
                                      low_counts_threshold, high_counts_threshold, 
                                      low_genes_threshold, high_genes_threshold, mt_threshold, file_path2 = file_path2)
cat("QC filtering completed successfully.\n")

# --- Step 5: Normalise and Scale the Data ---
cat("Step 5: Normalising and scaling the data\n")
sample_object <- run_normalisation(sample_object, file_path2 = file_path2)
cat("Normalisation and scaling completed successfully.\n")

# --- Step 6: Perform Cell Cycle Scoring ---
cat("Step 6: Cell cycle scoring\n")
sample_object <- cell_cycle_scoring(sample_object)
cat("Cell cycle scoring completed.\n")

# --- Step 7: Plot Variable Features ---
cat("Step 7: Plotting variable features\n")
sample_object <- generate_variable_features_plot(sample_object, sample_name, file_path2 = file_path2)
cat("Variable features plot generated successfully.\n")

# --- Step 8: Run Principal Component Analysis (PCA) ---
cat("Step 8: Running PCA\n")
sample_object <- run_pca(sample_object, sample_name, file_path2 = file_path2)
cat("PCA completed successfully.\n")

# Save the preprocessed object after PCA.
qsave(sample_object, file = file.path(sample_dir_data, paste0("1_", sample_name, "_preprocessed.qs")))
cat("Preprocessed data saved successfully.\n")

# --- Step 9: Determine the Number of PCs to Use ---
cat("Step 9: Determining number of PCs\n")
output1 <- determine_pcs(sample_object, sample_name, file_path2 = file_path2)
calc_pcs_variance <- output1$pcs_variance
calc_pcs_jackstraw <- output1$pcs_jackstraw
sample_object <- output1$sample_object
pcs <- max(calc_pcs_variance, calc_pcs_jackstraw)
cat(paste("Number of PCs determined:", pcs, "\n"))
qsave(sample_object, file = file.path(sample_dir_data, paste0("1_", sample_name, "_preprocessed.qs")))
cat("Preprocessed data saved successfully.\n")

# --- Step 10: Clustering & Clustree Plots ---
cat("Step 10: Performing clustering and generating clustree plots\n")
output2 <- perform_clustree(sample_object, pcs, sample_name)
plot1 <- output2$plot1
plot2 <- output2$plot2
sample_object <- output2$sample_object
# Run neighbor finding for the next clustering step.
sample_object <- FindNeighbors(sample_object, dims = 1:pcs)
cat("Clustering and clustree plots completed successfully.\n")

# --- Step 11: Apply Chosen Resolution for Clustering ---
cat("Step 11: Applying chosen resolution for clustering\n")
if (!exists("chosen_resolution")) {
  chosen_resolution <- as.numeric(readline(prompt = "Enter chosen resolution: "))
}
sample_object <- FindClusters(sample_object, reduction = "pca", dims = 1:pcs, resolution = chosen_resolution)
cat(paste("Clusters found with resolution:", chosen_resolution, "\n"))

# --- Step 12: Run Dimensionality Reduction (t-SNE and UMAP) ---
cat("Step 12: Performing dimensionality reduction (t-SNE and UMAP)\n")
sample_object <- run_dimred(sample_object, sample_name, pcs)
cat("t-SNE and UMAP completed successfully.\n")

# --- Step 13: Run Doublet Detection ---
cat("Step 13: Running DoubletFinder\n")
output3 <- doubletfinder(sample_object, sample_name, pcs, exp_multiplet)
sample_object <- output3$sample_object
pN <- output3$pN
pK <- output3$pK
doublets <- output3$doublets
singlets <- output3$singlets
cat("DoubletFinder completed successfully.\n")

# --- Step 14: Save the Preprocessed Seurat Object ---
cat("Step 14: Saving preprocessed data\n")
qsave(sample_object, file = file.path(sample_dir_data, paste0("2_", sample_name, "_preprocessed.qs")))
cat("Preprocessed data saved successfully.\n")

###########################################
# UPDATE METADATA WITH SAMPLE-SPECIFIC INFO
###########################################
metadata_file <- file.path(data_dir, "seurat_metadata.csv")
if (file.exists(metadata_file)) {
  seurat_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  if (sample_name %in% seurat_metadata$Sample_Name) {
    cat("Sample already exists in metadata. Updating the entry.\n")
    seurat_metadata <- seurat_metadata %>% filter(Sample_Name != sample_name)
  } else {
    cat("Sample does not exist in metadata. Adding new entry.\n")
  }
  new_entry <- data.frame(Sample_Name = sample_name, Sample_ID = sample_id, N_PCs_used = pcs, 
                          Chosen_Resolution = chosen_resolution, pN = pN, pK = pK, 
                          doublets = doublets, singlets = singlets, stringsAsFactors = FALSE)
  seurat_metadata <- rbind(seurat_metadata, new_entry)
  write.csv(seurat_metadata, file = metadata_file, row.names = FALSE)
  cat("Metadata updated successfully.\n")
} else {
  seurat_metadata <- data.frame(Sample_Name = character(), Sample_ID = character(), N_PCs_used = numeric(), 
                                Chosen_Resolution = numeric(), pN = numeric(), pK = numeric(), 
                                doublets = numeric(), singlets = numeric(), stringsAsFactors = FALSE)
  new_entry <- data.frame(Sample_Name = sample_name, Sample_ID = sample_id, N_PCs_used = pcs, 
                          Chosen_Resolution = chosen_resolution, pN = pN, pK = pK, 
                          doublets = doublets, singlets = singlets, stringsAsFactors = FALSE)
  seurat_metadata <- rbind(seurat_metadata, new_entry)
  write.csv(seurat_metadata, file = metadata_file, row.names = FALSE)
  cat("New metadata file created and updated successfully.\n")
}

###########################################
# GENE MARKERS IDENTIFICATION & HEATMAP GENERATION
###########################################
# Reload the latest preprocessed object.
sample_object <- qread(file.path(sample_dir_data, paste0("2_", sample_name, "_preprocessed.qs")))
# Identify marker genes for each cluster using Wilcoxon test and add percentage differences.
sample_id.markers <- FindAllMarkers(object = sample_object, verbose = TRUE, only.pos = TRUE, test.use = "wilcox") %>% Add_Pct_Diff()
sample_id.sig <- sample_id.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  ungroup()
# Save the significant marker genes to a CSV file 
write.csv(sample_id.sig, file = file.path(sample_dir_data, paste0(sample_id, "_significant_markers.csv")))

# For each cluster, keep only markers with avg_log2FC > 1 and take the top 20.
sample_id.top20 <- sample_id.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup()

# Save the top 20 marker genes for each cluster.
write.csv(sample_id.top20, file = file.path(sample_dir_data, paste0(sample_id, "_top20markergenes.csv")))

# Create and save a heatmap of the top marker genes.
plot1 <- DoHeatmap(sample_object, features = sample_id.top20$gene)
ggsave(filename = file.path(sample_dir_plots, paste0("13_", sample_name, "_MarkerGene_top20_HeatMap.pdf")),
       plot = plot1, width = 10, height = 20)

###########################################
# PLOTTING
###########################################
# Read in the Seurat object 
sample_object <- qread(file.path(sample_dir_data, paste0("2_", sample_name, "_preprocessed.qs")))

# UMAP plot by Seurat clusters.
sample_object <- SetIdent(sample_object, value = "seurat_clusters")
unique_clusters <- unique(sample_object@meta.data$seurat_clusters)
umap_plot_clusters <- dittoDimPlot(
  object = sample_object,
  var = "seurat_clusters",
  reduction.use = "umap",
  do.label = TRUE, labels.size = 4,
  main = "UMAP by Seurat Clusters"
)
print(umap_plot_clusters)
ggsave(filename = file.path(plots_dir, sample_name, paste0("14_", sample_name, "_BDGP_seurat_fixed.pdf")),
       plot = umap_plot_clusters, width = 8, height = 5, units = "in", dpi = 300)

# UMAP plot by timepoint.
sample_object <- SetIdent(sample_object, value = "timepoint")
umap_timepoint <- dittoDimPlot(
  object = sample_object,
  var = "timepoint",
  reduction.use = "umap",
  main = "UMAP by Timepoint"
)
print(umap_timepoint)
ggsave(filename = file.path(plots_dir, sample_name, paste0("15_", sample_name, "_Timepoint.pdf")),
       plot = umap_timepoint, width = 8, height = 5, units = "in", dpi = 300)

# UMAP plot by cell cycle phase.
sample_object <- SetIdent(sample_object, value = "Phase")
umap_phase <- dittoDimPlot(
  object = sample_object,
  var = "Phase",
  reduction.use = "umap",
  main = "UMAP by Cell Cycle Phase"
)
print(umap_phase)
ggsave(filename = file.path(plots_dir, sample_name, paste0("16_", sample_name, "_Phase.pdf")),
       plot = umap_phase, width = 8, height = 5, units = "in", dpi = 300)

cat("Annotation and plotting steps completed successfully.\n")
