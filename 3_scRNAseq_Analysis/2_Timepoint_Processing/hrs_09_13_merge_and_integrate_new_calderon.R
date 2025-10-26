# ====================================================
# Integrated Analysis Script – Drosophila scRNA-seq
# ====================================================

# ────────────────────────────────────────────────────
# 1. Environment setup
# ────────────────────────────────────────────────────
.libPaths('~/R/x86_64-pc-linux-gnu-library/4.5/')
options(future.globals.maxSize = 40000 * 1024^2)
setwd("/DataDrives/Drive2/Clifton/R_Projects/2025_Drosophila_scRNAseq_MonoamineSpecification/ANALYSIS/Step_2_Initial_preprocessing/Step_2.2_Timepoint_preprocessing/")

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
seurat_obj_1 <- qread(paste0("../Step_2.1_Sample_preprocessing/Data/Timepoint_Datasets/hrs_09_13_A6_Rep1/2_hrs_09_13_A6_Rep1_preprocessed.qs"))
seurat_obj_2 <- qread(paste0("../Step_2.1_Sample_preprocessing/Data/Timepoint_Datasets/hrs_09_13_A14_Rep2/2_hrs_09_13_A14_Rep2_preprocessed.qs"))

sample_name_1        <- "hrs_09_13_A6_Rep1"
sample_name_2        <- "hrs_09_13_A14_Rep2"
sample_name_main     <- "hrs_09_13"
sample_name          <- sample_name_main 

optimal_pcs               <- 95   # updated later by determine_pcs()
chosen_resolution         <- 2.9  # clustering before Harmony
chosen_resolution_harmony <- 2.3  # clustering after  Harmony

merge_plot_width      <- 26
merge_plot_height     <- 6
integrated_plot_width <- 26
integrated_plot_height<- 6

# Create project-specific output directories
data_dir <- file.path("Data/Timepoint_Datasets",  sample_name_main)
plot_dir <- file.path("Plots/Timepoint_Datasets", sample_name_main)
if (!dir.exists(data_dir))  dir.create(data_dir,  recursive = TRUE)
if (!dir.exists(plot_dir))  dir.create(plot_dir,  recursive = TRUE)

bdgp_file  <- file.path(file_path2, "Annotation_Tables/insitu_annot_hrs_09_13.csv")

# ====================================================
# Main Seurat analysis pipeline
# ====================================================

# ── Step 1. Merge Seurat objects ─────────────────────
print("Step 1: Merging Seurat objects")
merged_seurat  <- merge(x = seurat_obj_1, y = seurat_obj_2, project = sample_name_main)
merged_filename <- file.path(data_dir, paste0("1_", sample_name_main, "_merged.qs"))
qsave(merged_seurat, merged_filename)
print(paste("Merged Seurat object saved to:", merged_filename))

# ── Step 2. Normalise merged data ────────────────────
print("Step 2: Normalising merged data")
merged_seurat <- NormalizeData(merged_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# ── Step 3. Identify highly variable features ────────
print("Step 3: Identifying variable features")
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 5000)

# ── Step 4. Scale data ───────────────────────────────
print("Step 4: Scaling data")
all.genes   <- rownames(merged_seurat)
merged_seurat <- ScaleData(merged_seurat, features = all.genes)

# ── Step 5. Principal-component analysis ─────────────
print("Step 5: Running PCA")
merged_seurat <- RunPCA(merged_seurat, npcs = 200)
plot_elbow    <- ElbowPlot(merged_seurat, ndims = 200)
ggsave(file.path(plot_dir, paste0("1_", sample_name_main, "_ElbowPlot.png")),
       plot_elbow, width = 10, height = 6)
qsave(merged_seurat, merged_filename)

# ── Step 6. Determine the optimal number of PCs ──────
print("Step 6: Determining optimal PCs")
pcs_output    <- determine_pcs(merged_seurat, sample_name_main)
pcs_variance  <- pcs_output$pcs_variance
pcs_jackstraw <- pcs_output$pcs_jackstraw
optimal_pcs   <- max(pcs_variance, pcs_jackstraw)  # conservative choice
qsave(merged_seurat, merged_filename)

# ── Step 7. Clustering & clustree visualisation ──────
print("Step 7: Clustering and generating clustree plots")
output        <- perform_clustree(merged_seurat, optimal_pcs, sample_name_main)
merged_seurat <- output$sample_object
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:optimal_pcs)

# ── Step 8. Final clustering (pre-Harmony) ───────────
print("Step 8: Applying chosen resolution (pre-Harmony)")
merged_seurat <- FindClusters(merged_seurat, reduction = "pca",
                              dims = 1:optimal_pcs, resolution = chosen_resolution)

# ── Step 9. Dimensionality reduction (t-SNE/UMAP) ────
print("Step 9: Running t-SNE and UMAP")
merged_seurat <- run_dimred(merged_seurat, sample_name_main, optimal_pcs)
qsave(merged_seurat, merged_filename)

# ── Step 10. Batch integration with Harmony ──────────
print("Step 10: Running Harmony integration")
integrated_seurat <- RunHarmony(merged_seurat, group.by.vars = "batch",
                                reduction = "pca", dims.use = 1:optimal_pcs,
                                reduction.save = 'runharmony')
integrated_seurat <- JoinLayers(integrated_seurat)

# ── Step 11. Clustering after Harmony ────────────────
print("Step 11: Clustering post-Harmony")
output2          <- perform_clustree_harmony(integrated_seurat, optimal_pcs, sample_name_main)
integrated_seurat <- output2$sample_object
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:optimal_pcs)

# ── Step 12. Final clustering (post-Harmony) ─────────
print("Step 12: Applying chosen resolution (post-Harmony)")
integrated_seurat <- FindClusters(integrated_seurat, reduction = "pca",
                                  dims = 1:optimal_pcs, resolution = chosen_resolution_harmony)

# ── Step 13. Dimensionality reduction (post-Harmony) ─
print("Step 13: Running t-SNE and UMAP (post-Harmony)")
integrated_seurat <- run_dimred(integrated_seurat, sample_name_main, optimal_pcs)
qsave(integrated_seurat,
      file.path(data_dir, paste0("2_", sample_name_main, "_harmonyintegrated.qs")))

# ── Step 14. Update run-level metadata table ─────────
print("Step 14: Updating Seurat metadata table")
if (file.exists("Data/Timepoint_Datasets/seurat_metadata_merged.csv")) {
  seurat_metadata <- read.csv("Data/Timepoint_Datasets/seurat_metadata_merged.csv", stringsAsFactors = FALSE)
  if (sample_name_main %in% seurat_metadata$Sample_Name) {
    # Replace existing entry
    seurat_metadata <- seurat_metadata %>% filter(Sample_Name != sample_name_main)
  }
  new_entry <- data.frame(Sample_Name      = sample_name_main,
                          N_PCs_used       = optimal_pcs,
                          Chosen_Resolution= chosen_resolution_harmony,
                          stringsAsFactors = FALSE)
  seurat_metadata <- rbind(seurat_metadata, new_entry)
  write.csv(seurat_metadata, "Data/Timepoint_Datasets/seurat_metadata_merged.csv", row.names = FALSE)
  print("Metadata updated.")
} else {
  # Create new metadata file
  seurat_metadata <- data.frame(Sample_Name      = sample_name_main,
                                N_PCs_used       = optimal_pcs,
                                Chosen_Resolution= chosen_resolution_harmony,
                                stringsAsFactors = FALSE)
  write.csv(seurat_metadata, "Data/Timepoint_Datasets/seurat_metadata_merged.csv", row.names = FALSE)
  print("Metadata file created.")
}

# ── Step 15. Marker-gene analysis & heatmap ──────────
integrated_seurat <- qread(file.path(data_dir, paste0("2_", sample_name_main, "_harmonyintegrated.qs")))
sample_id.markers <- FindAllMarkers(integrated_seurat, verbose = TRUE, only.pos = TRUE, test.use = "wilcox") %>%
  Add_Pct_Diff()

# Keep significant markers (FDR < 0.05)
sample_id.markers.sig <- sample_id.markers[sample_id.markers$p_val_adj < 0.05, ]
write.csv(sample_id.markers.sig,
          file.path(data_dir, paste0(sample_name_main, "_markergenes_sig.csv")))

# Top-20 markers per cluster (log2FC > 1)
sample_id.top20 <- sample_id.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup()
write.csv(sample_id.top20,
          file.path(data_dir, paste0(sample_name_main, "_top20markergenes.csv")))

plot1 <- DoHeatmap(integrated_seurat, features = sample_id.top20$gene)
ggsave(file.path(plot_dir, paste0("13_", sample_name_main, "_MarkerGene_top20_HeatMap.png")),
       plot1, width = 10, height = 20)

qsave(integrated_seurat,
      file.path(data_dir, paste0("2_", sample_name_main, "_harmonyintegrated.qs")))

print("Integration and marker analysis completed.")

# ====================================================
# Annotation - BDGP enrichment analysis
# ====================================================

# (Optional) Fisher test for marker enrichment
# fisher_results <- perform_fishers_test(
#   sample_id            = sample_name_main,
#   seurat_markers_path  = file.path(data_dir, paste0(sample_name_main, "_top20markergenes.csv")),
#   custom_db_path       = "Annotation_Tables/insitu_annot_hrs_00_03.csv",
#   output_dir           = data_dir
# )

# The following code incorporates your new BDGP enrichment functions and steps.
# It annotates scRNA-seq clusters by testing for over-represented BDGP terms.
integrated_seurat <- qread(file = paste0(data_dir, "/2_", sample_name_main, "_harmonyintegrated.qs"))

# (A) Define the BDGP enrichment function
# BDGP enrichment functions-
# This functions runs a fisher test to find over-represented BDGP terms
# Input tables must have 2 columns: "cluster" and "term"
# For each cluster, a result table is saved in outdir

# (A) Define the BDGP enrichment function using a clean for-loop and lapply for Fisher tests.
enrichment.BDGP <- function(fg.table, bg.table, outdir) {
  unique_clusters <- sort(unique(fg.table$cluster))
  for (i in unique_clusters) {
    # Extract foreground terms for the current cluster.
    fg_terms <- fg.table$term[fg.table$cluster == i]
    bg_terms <- bg.table$term
    # Split comma-separated terms, trim whitespace, and get unique term list.
    terms_list <- strsplit(fg_terms, ",")
    all_terms <- unique(trimws(unlist(terms_list)))
    
    message(sprintf("Cluster %s: testing %d genes with %d unique terms. Background genes: %d",
                    i, length(fg_terms), length(all_terms), length(bg_terms)))
    
    # Run Fisher's exact test for each unique term.
    results <- lapply(all_terms, function(term) {
      fg_yes <- sum(grepl(term, fg_terms))
      fg_no <- length(fg_terms) - fg_yes
      bg_yes <- sum(grepl(term, bg_terms))
      bg_no <- length(bg_terms) - bg_yes
      
      # Build contingency matrix for the Fisher test.
      mat <- matrix(c(fg_yes, fg_no, bg_yes, bg_no), nrow = 2)
      ft <- fisher.test(mat, alternative = "greater")
      # Apply Bonferroni correction (multiply p-value by number of tests).
      corrected_p <- ft$p.value * length(all_terms)
      
      data.frame(term = term, corrected.pval = corrected_p, significant = corrected_p < 0.05,
                 stringsAsFactors = FALSE)
    })
    
    # Combine results and sort by corrected p-value.
    res <- do.call(rbind, results)
    res <- res[order(res$corrected.pval), ]
    # Save the enrichment results for the current cluster.
    write.table(res, file = file.path(outdir, paste0("cluster", i, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# (B) Prepare BDGP annotation table and input files.
tp <- sample_name_main  # Use the sample name as timepoint label, adjust if necessary.
# Read the BDGP annotation table; assume first column is gene names.
bdgp_table <- read.table(bdgp_file, quote = "", header = FALSE, sep = ",")
# Combine all annotation columns (except the gene name column) into one string.
temp <- bdgp_table[, colnames(bdgp_table) != 'V1', drop = FALSE]
temp[temp == ""] <- NA
library(tidyr)
bdgp_table$term <- unite(temp, 'term', sep = ',', na.rm = TRUE)$term

# (C) Create the background gene list from a previously processed object.
cell.data <- qread(paste0(data_dir, "/2_", sample_name_main, "_harmonyintegrated.qs"))
genes <- row.names(cell.data)
umi_count <- rowSums(cell.data)
genes_umi <- data.frame(genes = genes, umi_count = as.numeric(umi_count), row.names = genes)

# Filter genes with umi_count greater than 1.
bg <- genes_umi[genes_umi$umi_count > 1, ]
write.csv(bg, file = paste0("bg.", tp, ".csv"), row.names = FALSE)

# (D) Prepare the foreground gene list (marker genes).
library(data.table)
fg <- fread(file.path(paste0("Data/Timepoint_Datasets/", sample_name_main , "/", sample_name_main, "_markergenes_sig.csv")))
fg <- fg %>%
  dplyr::filter(p_val_adj < 0.05,
                avg_log2FC >= 0.25)

# Filter out unwanted (e.g., ribosomal) genes; ensure you have a file "ribo_genes.txt"
output_fg  <- paste0("fg.", tp, ".no.ub.txt")
output_bg  <- paste0("bg.", tp, ".no.ub.txt")
fg %>% filter(!grepl("^(RpL|RpS)", gene)) %>% write.table(output_fg, sep = "\t", quote = FALSE, row.names = FALSE)
bg %>% filter(!grepl("^(RpL|RpS)", genes)) %>% write.table(output_bg, sep = "\t", quote = FALSE, row.names = FALSE)

# (E) Remove fg genes from bg (to ensure non-DE background)
output_file  <- paste0("bg.", tp, ".no.ub.final.txt")
bg <- fread(output_bg)
fg <- fread(output_fg)
bg %>% filter(!(genes %in% fg$gene)) %>% write.table(output_file, sep = "\t", quote = FALSE, row.names = FALSE)
bg <- fread(output_file)
bg$cluster <- "BG"
bg <- bg[, c(1,3)]
bg %>% write.table(output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# (F) Assign BDGP terms to genes in fg and bg by matching with the BDGP annotation table.
bg_genes <- read.table(output_file, header = TRUE, stringsAsFactors = FALSE)
fg_genes <- read.table(output_fg, header = TRUE, stringsAsFactors = FALSE)
bg_genes$term <- bdgp_table$term[match(bg_genes$gene, bdgp_table$V1)]
fg_genes$term <- bdgp_table$term[match(fg_genes$gene, bdgp_table$V1)]

# Keep only rows with non-missing term information and required columns.
fg_final <- fg_genes[!is.na(fg_genes$term), c("cluster", "term")]
fg_final$term <- as.character(fg_final$term)
bg_final <- bg_genes[!is.na(fg_genes$term), c("cluster", "term")]
bg_final$term <- as.character(bg_final$term)

# (G) Create an output directory for BDGP enrichment results.
bdgp_outdir <- paste0(data_dir, "/", tp, "_BDGP_Enrichment")
if (!dir.exists(bdgp_outdir)) dir.create(bdgp_outdir, recursive = TRUE)

# (H) Run BDGP enrichment analysis by cluster.
enrichment.BDGP(fg_final, bg_final, bdgp_outdir)

# (I) Optionally, compile results from individual cluster files into a single CSV.
dat <- data.frame(row.names = 1:10)
cluster_files <- sort(list.files(bdgp_outdir, pattern = ".txt", full.names = TRUE))
for (fl in cluster_files) {
  dat0 <- read.delim(fl)
  if (nrow(dat0) < 10) {
    dat0[(nrow(dat0) + 1):(10), ] <- NA
  }
  colnames(dat0) <- c(fl, "corrected.pval", "significant")
  dat <- cbind(dat, dat0[1:10, 1:2])
}
write.csv(dat, file = paste0(data_dir, "/", tp, ".bdgp.enrich.final.csv"))

# (J) Optionally, save the top 20 foreground genes per cluster (with BDGP terms)
output_file_top20 <- paste0(data_dir, "/top20.genes.", tp, ".csv")
fg_genes %>% arrange(p_val_adj) %>% group_by(cluster) %>% slice(1:20) %>% write.csv(output_file_top20, quote = FALSE, row.names = FALSE)

# List the enrichment result files (assuming they are named like "cluster1.txt", "cluster2.txt", etc.)
cluster_files <- list.files(bdgp_outdir, pattern = "^cluster.*\\.txt$", full.names = TRUE)

# Read each file, add a cluster column extracted from the filename, and combine all data frames.
combined_enrichment <- do.call(rbind, lapply(cluster_files, function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Extract the cluster ID from the filename, e.g. "cluster3.txt" becomes "3"
  cluster_id <- sub("cluster(.*)\\.txt", "\\1", basename(file))
  df$cluster <- cluster_id
  df
}))

# Optionally, order the combined table by cluster and by corrected p-value.
combined_enrichment <- combined_enrichment[order(combined_enrichment$cluster, combined_enrichment$corrected.pval), ]

# Write the combined table to a CSV file.
output_combined <- paste0(data_dir, "/", tp, ".bdgp.enrich.combined.csv")
write.csv(combined_enrichment, file = output_combined, row.names = FALSE)

message("Combined enrichment table created and saved as: ", output_combined)
print("BDGP enrichment analysis completed successfully.")


# -----------------------------
# (Remaining Annotation and Plotting Steps)
# -----------------------------
# Define file paths
input_file <- file.path(data_dir, paste0("2_", sample_name_main, "_harmonyintegrated.qs"))
output_file <- file.path(data_dir, paste0("3_", sample_name_main, "_harmonyintegratedannotated.qs"))
color_csv <- "../Supplementary_Data/label2colour_mapping.csv"
plot_dir <- file.path("Plots/Timepoint_Datasets", sample_name_main)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load and Annotate Clusters ---
integrated_seurat <- qread(input_file)
integrated_seurat <- SetIdent(integrated_seurat, value = "seurat_clusters")
integrated_seurat@meta.data$celltype <- Idents(integrated_seurat)

# Define annotations
new_annotations <- c(
  "brain primordium",
  "brain primordium",
  "anterior midgut primordium",
  "amnioserosa",
  "embryonic dorsal epidermis",
  "anterior midgut primordium",
  "anterior midgut primordium",
  "anterior midgut primordium",
  "brain primordium",
  "anterior midgut primordium",
  "amnioserosa",
  "anterior midgut primordium",
  "anterior midgut primordium",
  "anterior midgut primordium",
  "anterior midgut primordium",
  "anterior midgut primordium",
  "amnioserosa",
  "brain primordium",
  "anterior midgut primordium",
  "brain primordium",
  "anterior midgut primordium",
  "anterior midgut primordium",
  "brain primordium",
  "anterior midgut primordium",
  "brain primordium",
  "anterior midgut primordium",
  "anterior midgut primordium",
  "brain primordium",
  "anterior midgut primordium",
  "anterior midgut primordium",
  "anterior midgut primordium",
  "brain primordium",
  "anterior midgut primordium",
  "brain primordium",
  "anterior midgut primordium",
  "anterior midgut primordium"
)

# Apply annotations
cluster_ids <- as.numeric(as.character(Idents(integrated_seurat)))
if (min(cluster_ids, na.rm = TRUE) == 0) cluster_ids <- cluster_ids + 1

stopifnot(length(new_annotations) >= max(cluster_ids))
integrated_seurat@meta.data$celltype <- new_annotations[cluster_ids]
integrated_seurat@meta.data$celltype_unique <- paste(
  integrated_seurat@meta.data$seurat_clusters,
  integrated_seurat@meta.data$celltype,
  sep = "_"
)

# --- Load Color Mapping ---
annotation_table <- read.csv(color_csv)
integrated_seurat <- SetIdent(integrated_seurat, value = "celltype")
valid_labels <- levels(Idents(integrated_seurat))
color_mapping <- setNames(
  annotation_table$Color[annotation_table$Label %in% valid_labels],
  annotation_table$Label[annotation_table$Label %in% valid_labels]
)
cell_colors <- color_mapping[as.character(integrated_seurat@meta.data$celltype)]
integrated_seurat@meta.data$celltype_color <- cell_colors

# --- Helper: Plot Function ---
save_umap_plot <- function(var, color_panel = NULL, file_suffix, legend = TRUE, width = 8, height = 5) {
  plot <- dittoDimPlot(
    object = integrated_seurat,
    var = var,
    reduction.use = "umap",
    color.panel = color_panel,
    main = paste("UMAP by", var)
  )
  if (!legend) plot <- plot + theme(legend.position = "none")
  ggsave(
    filename = file.path(plot_dir, paste0("14_", sample_name_main, "_", file_suffix, ".pdf")),
    plot = plot,
    width = width, height = height, units = "in", dpi = 300
  )
}

# --- Generate UMAPs ---
save_umap_plot("celltype", color_mapping, "BDGP_CellType")
save_umap_plot("celltype", color_mapping, "BDGP_CellTypenolegend", legend = FALSE, width = 6)

# Plot with unique celltype labels
integrated_seurat <- SetIdent(integrated_seurat, value = "celltype_unique")
unique_labels <- unique(integrated_seurat@meta.data$celltype_unique)
celltype_parts <- sub(".*_", "", unique_labels)
color_mapping_unique <- setNames(color_mapping[celltype_parts], unique_labels)
save_umap_plot("celltype_unique", color_mapping_unique, "BDGP_CellTypeUnique_Custom", width = 14)

# Plot by seurat_clusters with labels
integrated_seurat <- SetIdent(integrated_seurat, value = "seurat_clusters")
cluster_map <- setNames(
  integrated_seurat@meta.data$celltype[match(unique(integrated_seurat@meta.data$seurat_clusters), integrated_seurat@meta.data$seurat_clusters)],
  unique(integrated_seurat@meta.data$seurat_clusters)
)
cluster_colors <- setNames(color_mapping[cluster_map], names(cluster_map))

umap_plot_clusters <- dittoDimPlot(
  object = integrated_seurat,
  var = "seurat_clusters",
  reduction.use = "umap",
  color.panel = cluster_colors,
  do.label = TRUE, labels.size = 4,
  main = "UMAP by Seurat Clusters"
)
ggsave(
  filename = file.path(plot_dir, paste0("14_", sample_name_main, "_BDGP_seurat_fixed.pdf")),
  plot = umap_plot_clusters,
  width = 8, height = 5, units = "in", dpi = 300
)

# --- Plot by Timepoint, Phase, Batch ---
library(scales) 

for (id in c("Phase", "batch")) {
  integrated_seurat <- SetIdent(integrated_seurat, value = id)
  
  n_groups <- length(levels(integrated_seurat))  # how many colours needed
  pal      <- hue_pal()(n_groups)               # default ggplot palette
  
  save_umap_plot(id, pal, paste0(id))
}

# --- Load Colour Mapping ---
timepoint_color_csv <- "../Supplementary_Data/timepoint2colour_mapping.csv"
annotation_table2   <- read.csv(timepoint_color_csv)

integrated_seurat <- SetIdent(integrated_seurat, value = "timepoint")

valid_labels <- levels(Idents(integrated_seurat))   # current factor levels

color_mapping2 <- setNames(
  annotation_table2$Color[annotation_table2$Label %in% valid_labels],
  annotation_table2$Label[annotation_table2$Label %in% valid_labels]
)

# optional: warn if something is missing
missing <- setdiff(valid_labels, names(color_mapping2))
if (length(missing) > 0)
  warning("No colour specified for: ", paste(missing, collapse = ", "))

# call the helper with the CORRECT mapping
save_umap_plot("timepoint", color_mapping2, "Timepoint")


# --- Save Seurat Object ---
qsave(integrated_seurat, output_file)
print("Pipeline finished successfully.")
