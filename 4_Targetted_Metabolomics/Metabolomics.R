
#############################################################
# METABOLOMICS DATA ANALYSIS SCRIPT
# Author: Clifton Lewis
# Date: 18 Mar 2025
#
# Description:
# This script performs the following:
#   1. Loads and cleans the data, separating raw and log-transformed
#      datasets along with their respective metadata.
#   2. Performs PCA on both datasets and generates scree and score plots.
#   3. Runs statistical tests (ANOVA or Kruskal-Wallis) for each metabolite,
#      computes effect sizes, and produces annotated boxplots.
#   4. Creates supervised and unsupervised clustering heatmaps (using raw data).
#
# Interpretation Guidelines:
#   - PCA plots: Look at the scree plot to determine how many principal
#     components explain the variance. In the score plots, examine sample
#     clustering by group, batch, or injection order.
#   - Statistical tests: The printed p-values and effect sizes indicate which
#     metabolites differ significantly across timepoints (Group). ANOVA results
#     use partial eta squared, while Kruskal-Wallis results use a non-parametric
#     effect size.
#   - Heatmaps: Supervised heatmaps (with group annotations) help visualize
#     clustering patterns of metabolites. Unsupervised heatmaps show inherent
#     clustering without any external labels.
#############################################################

# -------------------------------
# 1. Load Libraries
# -------------------------------
library(tidyverse)
library(limma)
library(ggplot2)
library(pheatmap)
library(sva)
library(rstatix)
library(mixOmics)
library(splines)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(ggridges)
library(ggsci)
library(ggrepel)
library(factoextra)
library(corrplot)
library(dunn.test)
library(car)
library(ggpubr)
# -------------------------------
# 2. Data Input and Cleaning
# -------------------------------

# STEP 1: Load the Data
raw_data_curated <- read.delim("01_raw_data/Raw_Data/RawData_curated_R_w_QC.csv", 
                               header = TRUE, sep = ",", 
                               stringsAsFactors = FALSE, check.names = FALSE)

# STEP 2: Extract Metadata
# Assumption: Last 3 rows contain metadata.
metadata <- raw_data_curated[(nrow(raw_data_curated)-2):nrow(raw_data_curated), ]
# Remove metadata rows from main dataset:
data_main <- raw_data_curated[1:(nrow(raw_data_curated)-3), ]

# STEP 3: Transpose and Format Metadata
metadata <- as.data.frame(t(metadata))
colnames(metadata) <- metadata[1, ]  # First row becomes column names.
metadata <- metadata[-1, ]           # Remove redundant first row.
rownames(metadata) <- colnames(data_main)[-1]  # Use sample names (excluding first column)

# Convert selected metadata columns to factors
metadata$Injection_order <- as.factor(metadata$Injection_order)
metadata$Batch_ID <- as.factor(metadata$Batch_ID)

# STEP 4: Transpose the Main Data Matrix so that samples are rows.
data_main <- as.data.frame(t(data_main))
colnames(data_main) <- data_main[1, ]
data_main <- data_main[-1, ]

# STEP 5: Convert Data Matrix to Numeric
data_numeric <- as.data.frame(lapply(data_main, function(x) as.numeric(as.character(x))))
rownames(data_numeric) <- rownames(data_main)

# STEP 6: Create Combined Data with Both Raw and Log-Transformed Values
# Log-transform data using ln(x+1) to avoid log(0)
log_transformed <- log(data_numeric + 1)
# Append log-transformed columns with "_log" suffix
data_combined <- cbind(data_numeric, setNames(log_transformed, paste0(colnames(log_transformed), "_log")))

# STEP 7: Remove Specific QC Samples ("QC 1" and "QC 1C 2")
samples_to_remove <- grepl("QC 1$|QC 1C 2$", rownames(data_numeric), ignore.case = TRUE)
data_numeric <- data_numeric[!samples_to_remove, ]
metadata <- metadata[!samples_to_remove, ]

samples_to_remove <- grepl("QC 1$|QC 1C 2$", rownames(log_transformed), ignore.case = TRUE)
log_transformed <- log_transformed[!samples_to_remove, ]
metadata <- metadata[!samples_to_remove, ]

# STEP 8: Create Datasets with and without QC & Blank Samples
# Define sample filters based on sample names
samples_no_blank <- !grepl("Blank", rownames(data_numeric), ignore.case = TRUE)
samples_no_qc_blank <- !grepl("QC|Blank", rownames(data_numeric), ignore.case = TRUE)

# Dataset with QC but without Blank samples:
data_w_qc <- data_numeric[samples_no_blank, ]
metadata_w_qc <- metadata[samples_no_blank, ]

# Dataset without QC and Blank samples:
data_clean <- data_w_qc[samples_no_qc_blank, ]
metadata_clean <- metadata_w_qc[samples_no_qc_blank, ]

# STEP 9: Ensure Timepoint Information is Correct in Metadata
metadata_w_qc <- metadata_w_qc %>%
  tibble::rownames_to_column("Sample")
metadata_clean <- metadata_clean %>%
  tibble::rownames_to_column("Sample")
metadata_w_qc$Group <- metadata_w_qc$Timepoint
metadata_clean$Group <- metadata_clean$Timepoint
metadata_w_qc$Group <- as.factor(metadata_w_qc$Group)
metadata_clean$Group <- as.factor(metadata_clean$Group)

# STEP 10: Create Log-Transformed Datasets and Metadata
data_w_qc_log <- log_transformed[samples_no_blank, ]
metadata_w_qc_log <- metadata[samples_no_blank, ]
data_clean_log <- data_w_qc_log[samples_no_qc_blank, ]
metadata_clean_log <- metadata_w_qc_log[samples_no_qc_blank, ]
metadata_w_qc_log <- metadata_w_qc_log %>%
  tibble::rownames_to_column("Sample")
metadata_clean_log <- metadata_clean_log %>%
  tibble::rownames_to_column("Sample")
metadata_w_qc_log$Group <- metadata_w_qc_log$Timepoint
metadata_clean_log$Group <- metadata_clean_log$Timepoint
metadata_w_qc_log$Group <- as.factor(metadata_w_qc_log$Group)
metadata_clean_log$Group <- as.factor(metadata_clean_log$Group)

# -------------------------------
# 3. Define a Consistent Color Palette
# -------------------------------
# We choose a palette that will be used across all plots.
my_palette <- c("#00AFBB", "#E7B800", "#FC4E07", "#33CA7F", "#F26CA7", "#5E4AE3", "#A72608")
my_palette_faded <- adjustcolor(my_palette, alpha.f = 0.2)

# Optionally, force the Group factor to have a specific order (adjust levels as needed):
metadata_clean$Group <- factor(metadata_clean$Group, 
                               levels = c("0-3hpf", "3-6hpf", "6-9hpf", "9-12hpf", "12-15hpf", "15-18hpf", "18-22hpf"))
metadata_clean_log$Group <- factor(metadata_clean_log$Group, 
                                   levels = c("0-3hpf", "3-6hpf", "6-9hpf", "9-12hpf", "12-15hpf", "15-18hpf", "18-22hpf"))

# -------------------------------
# 4. PCA Analysis
# -------------------------------
# PCA is used here to visualize overall variation and clustering of samples.

# Perform PCA on raw data
pca_raw <- prcomp(scale(data_clean), center = TRUE, scale. = TRUE)
# Perform PCA on log-transformed data
pca_log <- prcomp(scale(data_clean_log), center = TRUE, scale. = TRUE)

# ----- Scree Plot for Raw Data -----
scree_data_raw <- data.frame(PC = seq_along(pca_raw$sdev),
                             Variance = (pca_raw$sdev^2) / sum(pca_raw$sdev^2))
scree_data_raw$Percentage <- round(scree_data_raw$Variance * 100, 2)
scree_plot_raw <- ggplot(scree_data_raw[1:8,], aes(x = factor(PC), y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(x = PC, y = Variance), colour = "red", size = 1, group = 1) +
  geom_point(size = 3, colour = "red") +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5, size = 4) +
  theme_minimal() +
  labs(title = "Scree Plot of PCA (Raw Data)",
       x = "Principal Component", 
       y = "Percentage Variance Explained")
# Add annotations explaining the components
scree_plot_raw_explained <- scree_plot_raw +
  annotate("text", x = 4, y = max(scree_data_raw$Variance[1:8]) * 0.9, 
           label = "Bars: Variance explained by each PC", size = 4, hjust = 0) +
  annotate("text", x = 4, y = max(scree_data_raw$Variance[1:8]) * 0.8, 
           label = "Line: Cumulative variance explained", size = 4, hjust = 0) +
  annotate("text", x = 4, y = max(scree_data_raw$Variance[1:8]) * 0.7, 
           label = "Points: Individual PC variance", size = 4, hjust = 0)
print(scree_plot_raw)
ggsave(filename = "Figures/PCA_raw_scree_plot.png", scree_plot_raw)
ggsave(filename = "Figures/PCA_raw_scree_plot_explained.png", scree_plot_raw_explained)

# ----- Scree Plot for Log-Transformed Data -----
scree_data_log <- data.frame(PC = seq_along(pca_log$sdev),
                             Variance = (pca_log$sdev^2) / sum(pca_log$sdev^2))
scree_data_log$Percentage <- round(scree_data_log$Variance * 100, 2)
scree_plot_log <- ggplot(scree_data_log[1:8,], aes(x = factor(PC), y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(x = PC, y = Variance), colour = "red", size = 1, group = 1) +
  geom_point(size = 3, colour = "red") +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5, size = 4) +
  theme_minimal() +
  labs(title = "Scree Plot of PCA (Log-Transformed Data)",
       x = "Principal Component", 
       y = "Percentage Variance Explained")
scree_plot_log_explained <- scree_plot_log +
  annotate("text", x = 4, y = max(scree_data_log$Variance[1:8]) * 0.9, 
           label = "Bars: Variance explained by each PC", size = 4, hjust = 0) +
  annotate("text", x = 4, y = max(scree_data_log$Variance[1:8]) * 0.8, 
           label = "Line: Cumulative variance explained", size = 4, hjust = 0) +
  annotate("text", x = 4, y = max(scree_data_log$Variance[1:8]) * 0.7, 
           label = "Points: Individual PC variance", size = 4, hjust = 0)
print(scree_plot_log)
ggsave(filename = "Figures/PCA_log_scree_plot.png", scree_plot_log)
ggsave(filename = "Figures/PCA_log_scree_plot_explained.png", scree_plot_log_explained)

# ----- PCA Score Plots -----
# Calculate percentage variance explained for PC1 and PC2
pve <- (pca_raw$sdev^2) / sum(pca_raw$sdev^2)
pc1_var <- round(pve[1] * 100, 1)
pc2_var <- round(pve[2] * 100, 1)

# PCA Plot (Raw Data) colored by Group and shaped by Batch with variance percentages in axis labels
pca_plot_batch_raw <- ggplot(
  data.frame(
    PC1 = pca_raw$x[,1],
    PC2 = pca_raw$x[,2],
    Group = metadata_clean$Group,
    Batch = metadata_clean$Batch_ID,
    Sample = metadata_clean$Sample
  ),
  aes(x = PC1, y = PC2, colour = Group, shape = Batch)
) +
  geom_point(size = 3) +
  #stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, level = 0.95) +
  geom_text(aes(label = Sample), hjust = 0.5, vjust = -1, size = 3) +
  scale_colour_manual(values = my_palette) +
  scale_fill_manual(values = my_palette_faded) +
  theme_bw(base_size = 16) +  # Increase the base font size
  labs(
    title = "PCA Plot (Raw Data) - Batch",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)")
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),  # Larger plot title
    axis.title = element_text(size = 16),                   # Larger axis titles
    axis.text = element_text(size = 14)                     # Larger axis text labels
  )

print(pca_plot_batch_raw)
ggsave(filename = "Figures/PCA_raw_batch_plot.png", pca_plot_batch_raw, width = 8, height = 6)

# -------------------------------
# 5. Statistical Analysis and Annotated Boxplots
# -------------------------------
# This function performs a statistical test for each metabolite,
# checks assumptions (normality, homogeneity of variances), and then applies:
#   - ANOVA (if parametric assumptions hold)
#   - Kruskal-Wallis (if not)
# Post-hoc tests are performed when needed.
# Results are saved as a CSV and individual logs per metabolite.

analyse_metabolites <- function(data, metadata) {
  results <- data.frame()
  log_text <- "Metabolomics Statistical Analysis Log\n\n"
  
  for(metabolite in colnames(data)) {
    # Create a temporary data frame for the metabolite and relevant metadata
    test_data <- data.frame(
      Value = data[, metabolite],
      Timepoint = metadata$Timepoint,
      Batch = metadata$Batch_ID
    )
    
    # Check Normality (Shapiro-Wilk) and Homogeneity of Variance (Levene's Test)
    shapiro_test <- shapiro.test(test_data$Value)
    levene_test <- leveneTest(Value ~ Timepoint, data = test_data)
    
    # Choose test based on assumption results
    if(shapiro_test$p.value > 0.05 && levene_test$`Pr(>F)`[1] > 0.05) {
      # Parametric assumptions met: Use ANOVA
      model <- aov(Value ~ Timepoint, data = test_data)
      anova_result <- summary(model)[[1]]
      p_value <- anova_result$`Pr(>F)`[1]
      test_used <- "ANOVA"
      
      # Post-hoc analysis if significant
      if(p_value < 0.05) {
        posthoc <- TukeyHSD(model)
        posthoc_summary <- capture.output(print(posthoc))
      } else {
        posthoc_summary <- "No significant post-hoc results"
      }
    } else {
      # Use non-parametric Kruskal-Wallis test
      kw_result <- kruskal.test(Value ~ Timepoint, data = test_data)
      p_value <- kw_result$p.value
      test_used <- "Kruskal-Wallis"
      
      if(p_value < 0.05) {
        posthoc <- dunn.test::dunn.test(test_data$Value,
                                        test_data$Timepoint,
                                        method = "bonferroni")
        posthoc_summary <- capture.output(print(posthoc))
      } else {
        posthoc_summary <- "No significant post-hoc results"
      }
    }
    
    # Store the results
    results <- rbind(results, data.frame(
      Metabolite = metabolite,
      Test_Used = test_used,
      P_Value = p_value,
      Significant = p_value < 0.05,
      Normality_P = shapiro_test$p.value,
      Variance_P = levene_test$`Pr(>F)`[1]
    ))
    
    # Create and save individual log for the metabolite
    metabolite_log <- paste0(
      "Metabolite: ", metabolite, "\n",
      "Test Used: ", test_used, "\n",
      "P-Value: ", round(p_value, 5), "\n",
      "Significant: ", ifelse(p_value < 0.05, "Yes", "No"), "\n",
      "Shapiro-Wilk Normality P: ", round(shapiro_test$p.value, 5), "\n",
      "Levene's Test P: ", round(levene_test$`Pr(>F)`[1], 5), "\n",
      "Post-hoc Test Summary:\n", paste(posthoc_summary, collapse = "\n"), "\n\n"
    )
    write(metabolite_log, file = paste0("Figures/", metabolite, "_Analysis_Log.txt"))
  }
  
  # Adjust for multiple testing using the Benjamini-Hochberg method
  results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
  
  # Save a full analysis log and CSV file with the results
  write(log_text, file = "Figures/Metabolite_Analysis_Log.txt")
  write.csv(results, "Figures/Metabolite_Analysis_Results.csv", row.names = FALSE)
  
  return(results)
}

# Run the analysis function on the raw data
results_raw <- analyse_metabolites(data = data_clean, metadata = metadata_clean)

results_log <- analyse_metabolites(data = data_clean_log, metadata = metadata_clean_log)

# ----- Annotated Boxplots for Each Metabolite -----
# For each metabolite, generate a boxplot with annotations that include:
#   - The test used (ANOVA or Kruskal-Wallis)
#   - The p-value
#   - The effect size (partial eta squared for ANOVA; non-parametric effect size for Kruskal-Wallis)
library(rstatix)
## Raw values
# Set your detection limit (adjust this value as needed)
detection_limit <- 100000  
for(met in results_raw$Metabolite) {
  
  # Prepare data for plotting:
  plot_data <- data.frame(
    Value = data_clean[[met]],
    Group = metadata_clean$Group
  )
  
  # Identify test used for this metabolite:
  test_used <- results_raw$Test_Used[results_raw$Metabolite == met]
  
  # Calculate effect size:
  if(test_used == "ANOVA"){
    # Use anova_test() from rstatix (partial eta squared is reported as 'pes')
    aov_res <- anova_test(data = plot_data, dv = Value, between = Group)
    effect_size <- as.numeric(aov_res$pes[1])
  } else {
    # For Kruskal-Wallis, use kruskal_effsize() from rstatix
    kw_eff <- kruskal_effsize(data = plot_data, formula = Value ~ Group)
    effect_size <- as.numeric(kw_eff$effsize)
  }
  
  # Get the adjusted p-value (Benjamini-Hochberg) and round it
  adj_p <- round(results_raw$Adjusted_P_Value[results_raw$Metabolite == met], 4)
  
  # Create a boxplot annotated with the test, adjusted p-value, effect size, and a horizontal dashed red line for the detection limit.
  p <- ggboxplot(plot_data, 
                 x = "Group", 
                 y = "Value", 
                 color = "Group", 
                 palette = my_palette,
                 add = "jitter", 
                 add.params = list(size = 1.5, alpha = 0.7)) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1) +
    # Add a horizontal dashed red line for the detection limit and map its linetype to include it in the legend
    geom_hline(aes(yintercept = detection_limit, linetype = "Detection limit"), 
               color = "red", size = 1) +
    scale_linetype_manual(name = "", values = c("Detection limit" = "dashed")) +
    labs(title = paste(met),
         subtitle = paste("Test:", test_used, 
                          "| Adjusted p-value (BH):", adj_p,
                          "| Effect Size:", round(effect_size, 3))) +
    scale_y_continuous(labels = scales::comma) +  # force full number formatting
    theme_minimal(base_size = 16) +  # Set a larger base font size
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)
    )
  
  print(p)
  ggsave(filename = paste0("Figures/", met, "_Boxplot.png"), plot = p, width = 10, height = 8, units = "in")
}


## Log-transformed values
# Set your detection limit for log-transformed data (adjust this value as needed)
# For example, if the raw detection limit is 100000 and you use ln(x+1), then:
detection_limit_log <- log(100000 + 1)

for(met in results_log$Metabolite) {
  
  # Prepare data for plotting using log-transformed data:
  plot_data_log <- data.frame(
    Value = data_clean_log[[met]],
    Group = metadata_clean_log$Group
  )
  
  # Identify test used for this metabolite:
  test_used <- results_log$Test_Used[results_log$Metabolite == met]
  
  # Calculate effect size on log-transformed data:
  if(test_used == "ANOVA"){
    # Use anova_test() from rstatix on the log-transformed data
    aov_res_log <- anova_test(data = plot_data_log, dv = Value, between = Group)
    effect_size_log <- as.numeric(aov_res_log$pes[1])
  } else {
    # Use kruskal_effsize() from rstatix on the log-transformed data
    kw_eff_log <- kruskal_effsize(data = plot_data_log, formula = Value ~ Group)
    effect_size_log <- as.numeric(kw_eff_log$effsize)
  }
  
  # Get the adjusted p-value (BH) and round it
  adj_p_log <- round(results_log$Adjusted_P_Value[results_log$Metabolite == met], 4)
  
  # Create a boxplot for the log-transformed data annotated with:
  # - The test used, adjusted p-value, and effect size
  # - A horizontal dashed red line representing the detection limit
  p_log <- ggboxplot(plot_data_log, 
                     x = "Group", 
                     y = "Value", 
                     color = "Group", 
                     palette = my_palette,
                     add = "jitter", 
                     add.params = list(size = 1.5, alpha = 0.7)) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1) +
    # Add a horizontal dashed red line for the detection limit
    geom_hline(aes(yintercept = detection_limit_log, linetype = "Detection limit"), 
               color = "red", size = 1) +
    scale_linetype_manual(name = "", values = c("Detection limit" = "dashed")) +
    labs(title = paste(met, "(Log-transformed)"),
         subtitle = paste("Test:", test_used, 
                          "| Adjusted p-value (BH):", adj_p_log,
                          "| Effect Size:", round(effect_size_log, 3))) +
    scale_y_continuous(labels = scales::comma) +  # force full number formatting
    theme_minimal(base_size = 16) +  # Set a larger base font size
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)
    )
  
  print(p_log)
  ggsave(filename = paste0("Figures/", met, "_Boxplot_Log.png"), plot = p_log, width = 10, height = 8, units = "in")
}

### ALternative ploting
for (met in results_log$Metabolite) {
  
  # Prepare data frame
  plot_data_log <- data.frame(
    Value = data_clean_log[[met]],
    Group = metadata_clean_log$Group
  )
  
  # Build dot‐only, mean‐line plot
  p_log <- ggplot(plot_data_log, aes(x = Group, y = Value, color = Group)) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", aes(group = 1),
                 color = "black", size = 1) +
    stat_summary(fun = mean, geom = "point",
                 color = "black", size = 3, shape = 18) +
    geom_hline(aes(yintercept = detection_limit_log, linetype = "Detection limit"),
               color = "red", size = 1) +
    scale_linetype_manual(name = "", values = c("Detection limit" = "dashed")) +
    labs(
      title = paste0(met, " (Log-transformed)"),
      x     = NULL,
      y     = "Log(Normalised Peak Area + 1)"
    ) +
    scale_y_continuous(labels = scales::comma) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(size = 12, face = "bold"),
      axis.title    = element_text(size = 12),
      axis.text     = element_text(size = 10),
      legend.text   = element_text(size = 10),
      legend.position = "none" 
    )
  
  print(p_log)
  ggsave(
    filename = paste0("Figures/", met, "_Dots_Log.png"),
    plot     = p_log,
    width    = 6, height = 4, units = "in"
  )
}

# -------------------------------
# 6. Combined Line Plots with Error Bars
# -------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# Assume 'data_clean_log' (log-transformed values) and 'metadata_clean_log' (with Sample and Group columns) are available.
# detection_limit_log is already defined as: 
# detection_limit_log <- log(100000 + 1)

### Plot 1: Neurotransmitters 
compounds1 <- c("Histamine", "Serotonin", "Octopamine", "Dopamine", "Tyramine")
# Define specific colours for each neurotransmitter using your pre-defined palette:
nt_colors <- c("Histamine" = my_palette[1], 
               "Serotonin" = my_palette[2], 
               "Octopamine" = my_palette[3], 
               "Dopamine" = my_palette[4], 
               "Tyramine" = my_palette[5])

# Prepare data by merging log-transformed values with metadata:
data_plot1 <- data_clean_log %>% 
  dplyr::mutate(Sample = rownames(data_clean_log)) %>%
  dplyr::select(Sample, all_of(compounds1))
data_plot1 <- left_join(data_plot1, metadata_clean_log, by = "Sample")

# Reshape data to long format:
data_plot1_long <- data_plot1 %>%
  pivot_longer(cols = all_of(compounds1), names_to = "Neurotransmitter", values_to = "Value")

# Summarize: compute mean and standard error (SE) for each timepoint and compound
summary_plot1 <- data_plot1_long %>%
  group_by(Group, Neurotransmitter) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            N = n(),
            SE = SD / sqrt(N),
            .groups = "drop")

# Create the line plot with error bars and the detection limit:
plot1 <- ggplot(summary_plot1, aes(x = Group, y = Mean, group = Neurotransmitter, colour = Neurotransmitter)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  geom_hline(yintercept = detection_limit_log, linetype = "dashed", color = "red", size = 1) +
  scale_colour_manual(values = nt_colors) +
  labs(title = "Monoamines",
       y = "ln (Peak Area)",
       x = "Timepoint") +
  theme_minimal(base_size = 16)
print(plot1)
ggsave(filename = "Figures/Neurotransmitters_LinePlot.png", plot = plot1, width = 10, height = 8, units = "in")


### Plot 2: Selected Metabolites (L-DOPA, tryosine, tryptophan)
compounds2 <- c("L.DOPA", "Tyrosine", "Tryptophan")

# Define colours for these metabolites (adjust as needed), and update the key for L-DOPA:
met_colors <- c("L.DOPA" = my_palette[6],
                "Tyrosine" = my_palette[7],
                "Tryptophan" = my_palette[1])

# Prepare data:
data_plot2 <- data_clean_log %>% 
  dplyr::mutate(Sample = rownames(data_clean_log)) %>%
  dplyr::select(Sample, all_of(compounds2))
data_plot2 <- left_join(data_plot2, metadata_clean_log, by = "Sample")

# Reshape data to long format:
data_plot2_long <- data_plot2 %>%
  pivot_longer(cols = all_of(compounds2), names_to = "Metabolite", values_to = "Value")

# Summarize by Group and Compound:
summary_plot2 <- data_plot2_long %>%
  group_by(Group, Metabolite) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            N = n(),
            SE = SD / sqrt(N),
            .groups = "drop")

# Create the second plot:
plot2 <- ggplot(summary_plot2, aes(x = Group, y = Mean, group = Metabolite, colour = Metabolite)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  geom_hline(yintercept = detection_limit_log, linetype = "dashed", color = "red", size = 1) +
  scale_colour_manual(values = met_colors) +
  labs(title = "Precursor Metabolites",
       y = "ln (Peak Area)",
       x = "Timepoint") +
  theme_minimal(base_size = 16)
print(plot2)
ggsave(filename = "Figures/MetabolitesPrecursors_LinePlot.png", plot = plot2, width = 10, height = 8, units = "in")


library(tidyverse)
library(ggplot2)

#----------------------------------------------------------
# RLA Analysis Using Log-Transformed Values (with QC data)
#----------------------------------------------------------

# For this analysis we use the log-transformed dataset with QC samples:
# data_w_qc_log and metadata_w_qc_log

# Remove duplicate columns (keeping the first occurrence)
metadata_w_qc_log <- metadata_w_qc_log[, !duplicated(names(metadata_w_qc_log))]

# Now add the rownames as a column, only if "Sample" is not already present:
if (!("Sample" %in% colnames(metadata_w_qc_log))) {
  metadata_w_qc_log <- metadata_w_qc_log %>% tibble::rownames_to_column("Sample")
} else {
  message("Column 'Sample' already exists after duplicate removal. No further action taken.")
}


### (A) RLA Across All Samples (Overall Median)
# Compute RLA: subtract the overall median (across all samples) for each metabolite.
rla_data_all <- data_w_qc_log - apply(data_w_qc_log, 2, median, na.rm = TRUE)

# Convert to long format for plotting.
rla_long_all <- rla_data_all %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "RLA")

# Merge with metadata to include Group info.
rla_long_all <- left_join(rla_long_all, metadata_w_qc_log, by = "Sample")

# Create a boxplot for RLA across all samples.
p_rla_all <- ggplot(rla_long_all, aes(x = Sample, y = RLA, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "RLA Plot - All Samples (Log Transformed, with QC)",
       y = "Relative Log Abundance (Overall Median Centered)",
       x = "Sample") +
  scale_fill_manual(values = c(my_palette, "grey")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p_rla_all)
ggsave("Figures/RLA_All_Samples_LogQC.png", p_rla_all, width = 10, height = 8, units = "in")

### (B) RLA Within Each Group (Group-Specific Median)
# Copy the log-transformed QC dataset to compute group-specific RLA.
rla_data_within <- data_w_qc_log

# Get unique groups from metadata
groups <- unique(metadata_w_qc_log$Group)

# For each group, subtract the median (computed within that group) for each metabolite.
for(grp in groups) {
  sample_grp <- metadata_w_qc_log$Sample[metadata_w_qc_log$Group == grp]
  if(length(sample_grp) > 1) {
    rla_data_within[sample_grp, ] <- rla_data_within[sample_grp, ] - 
      apply(rla_data_within[sample_grp, ], 2, median, na.rm = TRUE)
  }
}

# Convert the group-adjusted RLA data to long format.
rla_long_within <- rla_data_within %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "RLA")

# Merge with metadata.
rla_long_within <- left_join(rla_long_within, metadata_w_qc_log, by = "Sample")

# Create a boxplot for RLA computed within groups.
p_rla_within <- ggplot(rla_long_within, aes(x = Sample, y = RLA, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "RLA Plot - Within Groups (Log Transformed, with QC)",
       y = "Relative Log Abundance (Within Group Centered)",
       x = "Sample") +
  scale_fill_manual(values = c(my_palette, "grey")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p_rla_within)
ggsave("Figures/RLA_Within_Groups_LogQC.png", p_rla_within, width = 10, height = 8, units = "in")


# Order metadata by Group and Sample (using the log-transformed metadata)
ordered_metadata <- metadata_clean_log %>%
  arrange(Group, Sample)

# Reorder the log-transformed data rows to match the ordered metadata.
# Assumes rownames of data_clean_log correspond to sample names.
ordered_data_log <- data_clean_log[ordered_metadata$Sample, ]

# Prepare row annotation using only the Group column.
annotation_row_ordered <- ordered_metadata %>%
  dplyr::select(Sample, Group) %>%
  tibble::column_to_rownames("Sample")

# Define annotation colors for Group using your custom palette.
# Ensure that the number of colors matches the number of unique Group levels.
group_levels <- levels(metadata_clean_log$Group)
annotation_colors <- list(Group = setNames(my_palette[seq_along(group_levels)], group_levels))

# Create the supervised heatmap.
p_supervised <- pheatmap(
  ordered_data_log,
  annotation_row = annotation_row_ordered,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,  # keep the forced order
  cluster_cols = FALSE,  
  scale = "column",
  main = "Supervised Heatmap (Ordered by Group and Sample)",
  fontsize = 8
)


print(p_supervised)
# -------------------------------
# 8. Clustering Heatmaps (Raw Data Only)
# -------------------------------
# Supervised heatmap: Annotate each sample by its Group.
annotation_df <- data.frame(Group = metadata_clean$Group)
rownames(annotation_df) <- metadata_clean$Sample

pheatmap(
  scale(data_clean),
  annotation_row = annotation_df,
  color = colorRampPalette(my_palette)(50),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "Supervised Heatmap - Metabolite Clustering"
)
ggsave("Figures/Supervised_Heatmap.png", width = 10, height = 8, units = "in")

# Unsupervised heatmap: No sample annotation.
pheatmap(
  scale(data_clean),
  color = colorRampPalette(my_palette)(50),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "Unsupervised Heatmap - Metabolite Clustering"
)
ggsave("Figures/Unsupervised_Heatmap.png", width = 10, height = 8, units = "in")
