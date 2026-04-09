# =============================================================================
# Assignment 4 - scRNA-seq Analysis
# Script 06: Cell Composition Analysis
# =============================================================================

# --- Install packages (run once, then comment out) ---------------------------
# install.packages("compositions")
# install.packages("patchwork")

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Load checkpoint from script 05
seu <- readRDS("data/seurat_checkpoint_05.rds")

# =============================================================================
# PART 1: CELL TYPE PROPORTIONS ACROSS TIMEPOINTS (RM tissue)
# =============================================================================

# --- 1. Calculate proportions -------------------------------------------------
# Focus on RM tissue, exclude unassigned replicates
rm_cells <- seu@meta.data %>%
  filter(organ_custom == "RM",
         replicate    != "unassigned") %>%
  group_by(time, mouse_id, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(time, mouse_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Order timepoints correctly
rm_cells$time <- factor(rm_cells$time,
                        levels = c("Naive", "D02", "D05", "D08", "D14"))

# --- 2. Stacked bar plot: mean proportions per timepoint ----------------------
rm_mean_props <- rm_cells %>%
  group_by(time, cell_type) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

p_stacked <- ggplot(rm_mean_props,
                    aes(x = time, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell type composition over infection (RM tissue)",
       x     = "Timepoint",
       y     = "Mean proportion",
       fill  = "Cell type") +
  theme_classic() +
  theme(legend.position  = "right",
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(0.4, "cm"),
        axis.text.x      = element_text(size = 10))

ggsave("figures/06_composition_stacked_RM.png",
       p_stacked, width = 12, height = 7, dpi = 150, bg = "white")

# --- 3. Focus plot: immune cell proportions over time in RM -------------------
# Zoom in on just the immune cell types to see infection dynamics clearly
immune_types <- c("Macrophages", "Monocytes", "Neutrophils",
                  "Immature neutrophils", "Dendritic cells",
                  "NK cells", "T cells", "B cells")

rm_immune <- rm_cells %>%
  filter(cell_type %in% immune_types)

p_immune_time <- ggplot(rm_immune,
                        aes(x = time, y = proportion,
                            color = cell_type, group = cell_type)) +
  stat_summary(fun = mean, geom = "line",  linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
  labs(title  = "Immune cell proportions over infection (RM tissue)",
       x      = "Timepoint",
       y      = "Mean proportion ± SE",
       color  = "Cell type") +
  theme_classic() +
  theme(legend.position = "right")

ggsave("figures/06_immune_proportions_RM.png",
       p_immune_time, width = 10, height = 6, dpi = 150, bg = "white")

# =============================================================================
# PART 2: COMPOSITION ACROSS ALL THREE TISSUES
# =============================================================================

# --- 4. Proportions by tissue (all timepoints combined) ----------------------
tissue_props <- seu@meta.data %>%
  filter(replicate != "unassigned") %>%
  group_by(organ_custom, mouse_id, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(organ_custom, mouse_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

tissue_mean <- tissue_props %>%
  group_by(organ_custom, cell_type) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

p_tissue <- ggplot(tissue_mean,
                   aes(x = organ_custom, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell type composition by tissue",
       x     = "Tissue",
       y     = "Mean proportion",
       fill  = "Cell type") +
  theme_classic() +
  theme(legend.position = "right",
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"))

ggsave("figures/06_composition_by_tissue.png",
       p_tissue, width = 10, height = 7, dpi = 150, bg = "white")

# =============================================================================
# PART 3: STATISTICAL TEST — does macrophage proportion change over time in RM?
# =============================================================================
# Use a simple linear model on proportions per replicate.
mac_props <- rm_cells %>%
  filter(cell_type == "Macrophages") %>%
  mutate(time_num = case_when(time == "Naive" ~ 0,
                              time == "D02"   ~ 2,
                              time == "D05"   ~ 5,
                              time == "D08"   ~ 8,
                              time == "D14"   ~ 14))

# Linear model: proportion ~ timepoint
lm_mac <- lm(proportion ~ time_num, data = mac_props)
cat("Linear model: macrophage proportion ~ timepoint (RM)\n")
print(summary(lm_mac))

# Plot macrophage proportions per replicate over time
p_mac_prop <- ggplot(mac_props,
                     aes(x = time, y = proportion,
                         color = mouse_id, group = mouse_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               color = "black", linewidth = 1.2, linetype = "dashed") +
  labs(title  = "Macrophage proportion over infection (RM tissue)",
       x      = "Timepoint",
       y      = "Proportion of all RM cells",
       color  = "Mouse replicate") +
  theme_classic()

ggsave("figures/06_macrophage_proportion_RM.png",
       p_mac_prop, width = 8, height = 5, dpi = 150, bg = "white")

# --- 5. Save checkpoint -------------------------------------------------------
saveRDS(seu, "data/seurat_checkpoint_06.rds")

cat("\nScript 06 complete.\n")
