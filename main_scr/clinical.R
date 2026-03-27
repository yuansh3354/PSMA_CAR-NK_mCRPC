### ---------------
### Creator: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
### Date: 2024-02-27
### ---------------

# Part 1: Environment Initialization and Data Loading
# Clear workspace and reset garbage collection
rm(list = ls())
gc(reset = TRUE, verbose = FALSE)

# Set parallel computing memory limits and create output directory
options(future.globals.maxSize = 1000 * 1024^3)
script_start_time <- Sys.time()
# Load required libraries
library(readxl)
library(ggpattern)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(coin)
library(ggpubr)
library(eoffice)
library(patchwork)

# Part 2: Seurat Object Processing and Cell Filtering
# Read single-cell TCR integrated data
seurat_obj = qread(qs_tcr_file, nthreads = 512)

# Update sub-cluster voting results and remove NA cells
seurat_obj$subCluster_voting = ifelse(
  seurat_obj$majority_voting == 'B',
  'B', seurat_obj$subCluster_voting
)
cells_to_keep <- rownames(seurat_obj@meta.data[!is.na(seurat_obj$subCluster_voting), ])
seurat_obj <- seurat_obj[, cells_to_keep]

# Identify and remove incorrect cell types outside experimental design
seurat_obj$err = paste0(seurat_obj$sort, '_', seurat_obj$majority_voting)
remove_err <- c(
  "CARNK_Myeloid", "CARNK_CD8T", "CARNK_CD4T", "CARNK_B",
  "TC_Myeloid", "TC_NK", "TC_B", "CAR1_B"
)
keep_idx <- rownames(seurat_obj@meta.data[!(seurat_obj$err %in% remove_err), ])
seurat_obj <- seurat_obj[, keep_idx]

# Part 4: PSA Waterfall Plot Analysis (Figure 1a)
# Read maximum PSA decrease data
plt.df$`PSA Max` = plt.df$`PSA Max` * 100

# Generate Waterfall Plot
p1a <- ggplot(plot_data, aes(x = PatientID, y = PSA_Change)) +
  geom_col(aes(fill = Dose_Factor), color = "black", width = 0.8) +
  facet_grid(. ~ LC_Label, scales = "free_x", space = "free_x") +
  geom_text(aes(label = sprintf("%.1f", PSA_Change), hjust = ifelse(PSA_Change >= 0, -0.2, 1.2)), 
            angle = 90, size = 3) +
  geom_hline(yintercept = 0, linewidth = 0.8) +
  geom_hline(yintercept = c(-30, -50), linetype = "dashed", color = c("gray50", "darkred")) +
  theme_classic() +
  labs(y = "Max PSA Decrease from Baseline (%)", x = NULL) +
  scale_y_continuous(limits = c(-100, 100), breaks = c(-100, -50, -30, 0, 50, 100)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank())

topptx(p1a, file.path(dir_name, "Figure1a.pptx"), width = 8, height = 4)

# Part 5: Clinical Swimmer Plot (Figure 1b)
# Process clinical event timeline
event_colors <- c(
  "Screen" = "#D3D3D3", "LC" = "#A9A9A9", "CAR NK CELL INFUSION" = "#000000",
  "Stable disease" = "#56B4E9", "Partial relief" = "#009E73", 
  "PSA progression" = "#D55E00", "Survival after progression" = "#CC79A7"
)

# Clean time data
df_processed <- plt.df_swimmer %>%
  mutate(TimePoint = case_when(
    Days == "D0" ~ 0,
    grepl("-D-", Days) ~ -1 * as.numeric(sub(".*-D-", "", Days)),
    TRUE ~ as.numeric(sub("D", "", Days))
  ))

# Construct segment data for visualization
df_post <- df_processed %>%
  filter(!grepl("-", Days)) %>%
  arrange(PatientID, TimePoint) %>%
  group_by(PatientID) %>%
  mutate(Start = TimePoint, End = lead(TimePoint), Next_Event = lead(Event),
         Plot_Event = ifelse(Days == "D0", Next_Event, Event)) %>%
  filter(!is.na(End))

p1b <- ggplot() +
  geom_segment(data = df_post, aes(x = Start, xend = End, y = PatientID, yend = PatientID, color = Plot_Event), size = 6) +
  geom_point(data = df_processed %>% filter(Days == "D0"), aes(x = 0, y = PatientID), shape = 18, size = 4) +
  scale_color_manual(values = event_colors, name = "Clinical Status") +
  theme_classic() + labs(title = "Swimmer Plot", x = "Days relative to CAR-NK Infusion", y = NULL)

topptx(p1b, file.path(dir_name, "Figure1b.pptx"), width = 8, height = 5)

# Generate Adverse Event comparison plot
p1c <- ggplot(df_safety_stats) +
  geom_bar(aes(x = as.numeric(factor(Event)) - 0.2, y = LD_All, fill = "LD"), stat = "identity", width = 0.35) +
  geom_bar(aes(x = as.numeric(factor(Event)) + 0.2, y = NLD_All, fill = "Non-LD"), stat = "identity", width = 0.35) +
  geom_text(aes(x = as.numeric(factor(Event)), y = pmax(LD_All, NLD_All) + 0.5, label = P_Label), size = 3.5, fontface = "bold") +
  scale_x_continuous(breaks = 1:nrow(df_safety_stats), labels = df_safety_stats$Event) +
  scale_fill_manual(values = c("LD" = "#1F78B4", "Non-LD" = "#E31A1C"), name = "Group") +
  labs(title = "Comparison of Adverse Events", x = "", y = "Number of Patients") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

topptx(p1c, file.path(dir_name, "Figure1c.pptx"), width = 8, height = 4)

# Part 7: PSA and Cell Expansion Trend Analysis

# Define trend plotting function (SR/GR/Total)
plot_expansion_trend <- function(data, y_var, y_label) {
  ggplot(data, aes(x = Day_Numeric, y = !!sym(y_var))) +
    geom_line(aes(group = PatientID, color = Group), linewidth = 0.5, alpha = 0.15) +
    stat_summary(fun = mean, geom = "line", aes(group = 1, color = "Total"), linewidth = 1.5, linetype = "dashed") +
    stat_summary(fun = mean, geom = "line", aes(group = Group, color = Group), linewidth = 1.25) +
    scale_color_manual(values = my_colors, name = "Response Group") +
    scale_x_continuous(breaks = c(0, 3, 7, 14, 21), labels = c("D0", "D3", "D7", "D14", "D21"), name = "Time (Days)") +
    labs(y = y_label) + theme_classic()
}

p_nk_expansion <- plot_expansion_trend(flow_df, "CAR/NK%", "CAR-NK Expansion (%)")
topptx(p_nk_expansion, file.path(dir_name, "CAR_NK_Expansion.pptx"), width = 5, height = 4)

# Execution finished, return final status
message("All English translated figures and analysis completed.")