# SPATIAL ANALYSIS — total mortality (stacked bar) + per-pixel quantile boxplot
# by region x landcover (Wetland / Non-wetland) x elevation (< 5 m / > 5 m)
# Author: Henry CH Yeung, UVA

library(arrow)
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(colorspace)

set.seed(123)

base_dir <- "/Volumes/Henry2/ghostForest_us/_publish_product"   # EDIT to your local path
setwd(base_dir)

CONF_THRESH_STR <- "conf0d5"   # selects parquet version (and tags outputs)

pq_dir  <- "_driver_product/_full"
fig_dir <- paste0("_output/1_spatialAnalysis_", CONF_THRESH_STR)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ---- Load & filter ----
pq_files <- list.files(pq_dir,
                        pattern = paste0("_", CONF_THRESH_STR, "_.*\\.parquet$"),
                        full.names = TRUE)
cat(paste0("Parquet files (", CONF_THRESH_STR, "): ", length(pq_files), "\n"))

df <- rbindlist(lapply(pq_files, function(f) {
  as.data.table(read_parquet(f, col_select = c(
    "region", "dem_3dep", "slope", "nlcd", "mean_mort_rate", "all_mort",
    "n_obs_det", "insect_mort"
  )))
}))

df <- df[!is.na(mean_mort_rate) & insect_mort %in% c(0, NA)]
df[, group := fcase(nlcd == 90,              "Wetland",
                    nlcd %in% c(41, 42, 43), "Non-wetland")]
df <- df[!is.na(group) & dem_3dep >= 0]
df[, elev_cat := fifelse(dem_3dep > 5, "> 5 m", "< 5 m")]
df[, group4   := factor(paste0(group, ": ", elev_cat),
                        levels = c("Wetland: < 5 m", "Wetland: > 5 m",
                                   "Non-wetland: < 5 m", "Non-wetland: > 5 m"))]
df[, region   := factor(fcase(tolower(region) == "pacific",  "Pacific",
                               tolower(region) == "atlantic", "Atlantic",
                               tolower(region) == "gulf",     "Gulf",
                               tolower(region) == "lakes",    "Lakes"),
                        levels = c("Lakes", "Pacific", "Gulf", "Atlantic"))]
df <- df[!is.na(region)]
df_base <- df
df      <- df[n_obs_det > 10 | mean_mort_rate == 0]

# ---- Palette: blues = Wetland, browns = Non-wetland, dark = > 5 m ----
colors <- c("Non-wetland: > 5 m" = "#d9b98a",
            "Non-wetland: < 5 m" = "#ba8a50",
            "Wetland: > 5 m"     = "#8ec7e4",
            "Wetland: < 5 m"     = "#2a9cd5")

# ---- Panel b: stacked bars, Wetland left / Non-wetland right per region ----
bar_df <- df[, .(mort_m         = sum(all_mort, na.rm = TRUE) / 1e6,
                 area_km2       = .N / 100,                                  # 1 pixel = 1 ha = 0.01 km^2
                 hotspot_km2    = sum(mean_mort_rate >= 5, na.rm = TRUE) / 100,
                 hotspot_mort_m = sum(all_mort[mean_mort_rate >= 5], na.rm = TRUE) / 1e6),
               by = .(region, group4, group)]
bar_df[, xpos := as.integer(region) + fifelse(group == "Wetland", -0.15, 0.15)]

pb <- ggplot(bar_df, aes(x = xpos, y = mort_m, fill = group4, color = group4)) +
  geom_col(position = "stack", width = 0.26, linewidth = 0.6, alpha = 0.7) +
  scale_fill_manual(values = colors, name = NULL) +
  scale_color_manual(values = colorspace::darken(colors, 0.4), guide = "none") +
  scale_x_continuous(breaks = seq_along(levels(df$region)), labels = levels(df$region)) +
  scale_y_continuous(breaks = pretty_breaks(n = 4)) +
  labs(x = NULL, y = "Mortality (million)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.text.x     = element_text(angle = 0, hjust = 0.5))

# ---- Panel c: quantile boxplot of mean_mort_rate per pixel ----
box_df <- df[, .(
  mean   = mean(mean_mort_rate,           na.rm = TRUE),
  ymin   = quantile(mean_mort_rate, 0.10, na.rm = TRUE),
  lower  = quantile(mean_mort_rate, 0.25, na.rm = TRUE),
  middle = quantile(mean_mort_rate, 0.50, na.rm = TRUE),
  upper  = quantile(mean_mort_rate, 0.75, na.rm = TRUE),
  ymax   = quantile(mean_mort_rate, 0.90, na.rm = TRUE)
), by = .(region, group4)]

pc <- ggplot(box_df, aes(x = region, fill = group4, color = group4,
                          group = interaction(region, group4))) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle,
                   upper = upper, ymax = ymax),
               stat = "identity",
               position = position_dodge(0.75), width = 0.6,
               linewidth = 0.6, alpha = 0.7, show.legend = FALSE) +
  geom_point(aes(y = mean, fill = group4), position = position_dodge(0.75),
             size = 1.5, shape = 21, stroke = 0.8) +
  scale_fill_manual( values = colors, guide = "none") +
  scale_color_manual(values = colorspace::darken(colors, 0.4), guide = "none") +
  coord_cartesian(ylim = c(0, max(box_df$ymax))) +
  scale_y_continuous(breaks = pretty_breaks(n = 4)) +
  labs(x = NULL, y = expression("Mortality rate (ha"^-1*"yr"^-1*")"), shape = NULL) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# ---- Supplementary table: forest area, mortality & hotspot stats, by
#      region x landcover x elevation ----
lc_elev_table <- df[, .(
  `Forest area (km2)`                = .N / 100,
  `Mortality (million)`              = sum(all_mort, na.rm = TRUE) / 1e6,
  `Mortality rate (trees ha-1 yr-1)` = mean(mean_mort_rate, na.rm = TRUE),
  `Hotspot area (km2)`               = sum(mean_mort_rate >= 5, na.rm = TRUE) / 100,
  `Mortality in hotspot (million)`   = sum(all_mort[mean_mort_rate >= 5], na.rm = TRUE) / 1e6
), by = .(region, group, elev_cat)]

# Pooled across all regions ("All" row block)
lc_elev_table_pooled <- df[, .(
  region                              = "All",
  `Forest area (km2)`                = .N / 100,
  `Mortality (million)`              = sum(all_mort, na.rm = TRUE) / 1e6,
  `Mortality rate (trees ha-1 yr-1)` = mean(mean_mort_rate, na.rm = TRUE),
  `Hotspot area (km2)`               = sum(mean_mort_rate >= 5, na.rm = TRUE) / 100,
  `Mortality in hotspot (million)`   = sum(all_mort[mean_mort_rate >= 5], na.rm = TRUE) / 1e6
), by = .(group, elev_cat)]

lc_elev_table <- rbind(lc_elev_table,
                       lc_elev_table_pooled[, .(region, group, elev_cat,
                                                 `Forest area (km2)`, `Mortality (million)`,
                                                 `Mortality rate (trees ha-1 yr-1)`, `Hotspot area (km2)`,
                                                 `Mortality in hotspot (million)`)])
lc_elev_table[, `Forest area in hotspot (%)` := `Hotspot area (km2)` / `Forest area (km2)` * 100]

setnames(lc_elev_table, c("region", "group", "elev_cat"), c("Region", "Landcover", "Elevation"))
lc_elev_table[, Region    := factor(Region,    levels = c("Pacific", "Gulf", "Atlantic", "Lakes", "All"))]
lc_elev_table[, Landcover := factor(Landcover, levels = c("Wetland", "Non-wetland"))]
lc_elev_table[, Elevation := factor(Elevation, levels = c("> 5 m", "< 5 m"))]
setorder(lc_elev_table, Region, Landcover, Elevation)
lc_elev_table[, `:=`(Region = as.character(Region), Landcover = as.character(Landcover),
                     Elevation = as.character(Elevation))]
setcolorder(lc_elev_table, c("Region", "Landcover", "Elevation", "Forest area (km2)",
                             "Mortality (million)", "Mortality rate (trees ha-1 yr-1)",
                             "Hotspot area (km2)", "Forest area in hotspot (%)",
                             "Mortality in hotspot (million)"))
print(lc_elev_table)
fwrite(lc_elev_table, file.path(fig_dir, "table_landcover_elevation_summary.csv"))

# ---- Summary stats: mean, median, SD of mortality rate, Wetland vs Non-wetland ----
stats_by_region <- df[, .(
  mean = mean(mean_mort_rate, na.rm = TRUE), median = median(mean_mort_rate, na.rm = TRUE),
  sd = sd(mean_mort_rate, na.rm = TRUE), n = .N
), by = .(region, group)]

stats_pooled <- df[, .(
  region = "All",
  mean = mean(mean_mort_rate, na.rm = TRUE), median = median(mean_mort_rate, na.rm = TRUE),
  sd = sd(mean_mort_rate, na.rm = TRUE), n = .N
), by = .(group)]

# Region only (Wetland + Non-wetland pooled)
stats_region_only <- df[, .(
  group = "All",
  mean = mean(mean_mort_rate, na.rm = TRUE), median = median(mean_mort_rate, na.rm = TRUE),
  sd = sd(mean_mort_rate, na.rm = TRUE), n = .N
), by = .(region)]

# Grand total (all regions, both groups pooled)
stats_grand_total <- df[, .(
  region = "All", group = "All",
  mean = mean(mean_mort_rate, na.rm = TRUE), median = median(mean_mort_rate, na.rm = TRUE),
  sd = sd(mean_mort_rate, na.rm = TRUE), n = .N
)]

mort_summary <- rbind(stats_by_region,
                      stats_pooled[,      .(region, group, mean, median, sd, n)],
                      stats_region_only[, .(region, group, mean, median, sd, n)],
                      stats_grand_total)
mort_summary[, se := sd / sqrt(n)]
setcolorder(mort_summary, c("region", "group", "mean", "median", "sd", "se", "n"))
print(mort_summary)
fwrite(mort_summary, file.path(fig_dir, "stats_mortality_summary.csv"))

# ---- Pacific vs Non-Pacific summary stats (same mean/median/sd/se/n pattern) ----
df[, pacific_grp := fifelse(region == "Pacific", "Pacific", "Non-Pacific")]

stats_pacific_by_group <- df[, .(
  mean = mean(mean_mort_rate, na.rm = TRUE), median = median(mean_mort_rate, na.rm = TRUE),
  sd = sd(mean_mort_rate, na.rm = TRUE), n = .N
), by = .(pacific_grp, group)]

stats_pacific_pooled <- df[, .(
  group = "All",
  mean = mean(mean_mort_rate, na.rm = TRUE), median = median(mean_mort_rate, na.rm = TRUE),
  sd = sd(mean_mort_rate, na.rm = TRUE), n = .N
), by = .(pacific_grp)]

mort_summary_pacific <- rbind(stats_pacific_by_group,
                              stats_pacific_pooled[, .(pacific_grp, group, mean, median, sd, n)])
mort_summary_pacific[, se := sd / sqrt(n)]
setnames(mort_summary_pacific, "pacific_grp", "region")
setcolorder(mort_summary_pacific, c("region", "group", "mean", "median", "sd", "se", "n"))
print(mort_summary_pacific)
fwrite(mort_summary_pacific, file.path(fig_dir, "stats_mortality_summary_pacific.csv"))

# ---- Export bar/box data ----
fwrite(bar_df[, .(region, group4, mort_m, area_km2, hotspot_km2, hotspot_mort_m)],
       file.path(fig_dir, "stats_bar.csv"))
fwrite(box_df, file.path(fig_dir, "stats_box.csv"))

# ---- Legend (separate) ----
legend_plot <- ggplot(bar_df, aes(x = xpos, y = mort_m, fill = group4, color = group4)) +
  geom_col() +
  scale_fill_manual(values = colors, name = NULL) +
  scale_color_manual(values = colorspace::darken(colors, 0.4), guide = "none") +
  theme_void() +
  theme(legend.position       = "right",
        legend.key.size       = unit(10, "pt"),
        legend.text           = element_text(size = 11),
        legend.key.spacing.y  = unit(4, "pt"))

legend_grob <- cowplot::get_legend(legend_plot)
ggsave(file.path(fig_dir, "legend.pdf"), legend_grob, width = 1.75, height = 1.75, dpi = 300)
ggsave(file.path(fig_dir, "legend.png"), legend_grob, width = 1.75, height = 1.75, dpi = 300)

# ---- Combine panels b + c & save ----
g <- pb + pc + plot_layout(widths = c(1, 1.5))
ggsave(file.path(fig_dir, "spatial_mort_bc.pdf"), g, width = 12, height = 3.5, dpi = 300)
ggsave(file.path(fig_dir, "spatial_mort_bc.png"), g, width = 12, height = 3.5, dpi = 300)
plot(g)
cat("Saved to:", fig_dir, "\n")

# ================================================================
# Histograms of mortality distributions
# ================================================================
fill_col <- "#e96a80"
line_col <- colorspace::darken(fill_col, 0.4)
group_colors <- c("Non-wetland" = "#ba8a50", "Wetland" = "#2a9cd5")

# ---- Cumulative mortality (all_mort): overall ----
x_label_mort   <- expression(Cumulative ~ mortality ~ (ha^{-1}))
hotspot_cum_thresh <- 60   # 5 trees ha-1 yr-1, expressed as cumulative mortality over the study period

all_pct_overall <- mean(df$all_mort <= hotspot_cum_thresh, na.rm = TRUE) * 100
all_pct_label   <- paste0(round(all_pct_overall), "th percentile")


# ---- Cumulative mortality: Wetland vs Non-wetland ----
g_hist_all_group <- ggplot(df[all_mort > 0], aes(x = all_mort, fill = group, color = group)) +
  geom_histogram(alpha = 0.25, bins = 30, position = "identity") +
  geom_vline(xintercept = hotspot_cum_thresh, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  annotate("text", x = hotspot_cum_thresh, y = Inf, label = all_pct_label,
           hjust = 1.1, vjust = 1.2, size = 5, angle = 90, color = "grey30") +
  scale_x_sqrt(breaks = c(0, 10, 50, 100, 250, 500)) +
  scale_y_continuous(labels = function(x) x / 100) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_color_manual(values = colorspace::darken(group_colors, 0.4), guide = "none") +
  labs(x = x_label_mort, y = expression("Forest area (km"^2*")")) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.85, 0.85))
ggsave(file.path(fig_dir, "hist_allMort_overall_byGroup.pdf"), g_hist_all_group, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "hist_allMort_overall_byGroup.png"), g_hist_all_group, width = 6, height = 5, dpi = 300)
plot(g_hist_all_group)

# ---- Cumulative mortality: Wetland vs Non-wetland, faceted by region ----
all_pct_by_region <- df[all_mort > 0, .(pct = mean(all_mort <= hotspot_cum_thresh, na.rm = TRUE) * 100),
                        by = region]
all_pct_by_region[, label := paste0(round(pct), "th percentile")]

g_hist_all_group_region <- ggplot(df[all_mort > 0], aes(x = all_mort, fill = group, color = group)) +
  geom_histogram(alpha = 0.25, bins = 30, position = "identity") +
  geom_vline(xintercept = hotspot_cum_thresh, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_text(data = all_pct_by_region, inherit.aes = FALSE,
            aes(x = hotspot_cum_thresh, y = Inf, label = label),
            hjust = 1.1, vjust = 1.2, size = 4, angle = 90, color = "grey30") +
  scale_x_sqrt(breaks = c(0, 10, 50, 100, 250, 500)) +
  scale_y_continuous(labels = function(x) x / 100) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_color_manual(values = colorspace::darken(group_colors, 0.4), guide = "none") +
  facet_wrap(~ region, nrow = 1) +
  labs(x = x_label_mort, y = expression("Forest area (km"^2*")")) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(), legend.position = "bottom")
ggsave(file.path(fig_dir, "hist_allMort_overall_byGroup_region.pdf"),
       g_hist_all_group_region, width = 14, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "hist_allMort_overall_byGroup_region.png"),
       g_hist_all_group_region, width = 14, height = 4.5, dpi = 300)
plot(g_hist_all_group_region)

# ---- Mortality rate: overall ----
x_label_rate <- expression(Mortality ~ rate ~ (ha^{-1} * yr^{-1}))
hotspot_rate_thresh <- 5

rate_pct_overall <- mean(df$mean_mort_rate <= hotspot_rate_thresh, na.rm = TRUE) * 100
rate_pct_label   <- paste0(round(rate_pct_overall), "th percentile")

# ---- Mortality rate: Wetland vs Non-wetland ----
g_hist_rate_group <- ggplot(df[mean_mort_rate > 0], aes(x = mean_mort_rate, fill = group, color = group)) +
  geom_histogram(alpha = 0.25, bins = 30, position = "identity") +
  geom_vline(xintercept = hotspot_rate_thresh, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  annotate("text", x = hotspot_rate_thresh, y = Inf, label = rate_pct_label,
           hjust = 1.1, vjust = 1.2, size = 5, angle = 90, color = "grey30") +
  scale_x_sqrt(breaks = c(0, 1, 5, 10, 25, 50)) +
  scale_y_continuous(labels = function(x) x / 100) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_color_manual(values = colorspace::darken(group_colors, 0.4), guide = "none") +
  labs(x = x_label_rate, y = expression("Forest area (km"^2*")")) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.85, 0.85))
ggsave(file.path(fig_dir, "hist_meanMortRate_overall_byGroup.pdf"), g_hist_rate_group, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "hist_meanMortRate_overall_byGroup.png"), g_hist_rate_group, width = 6, height = 5, dpi = 300)
plot(g_hist_rate_group)

# ---- Mortality rate: Wetland vs Non-wetland, faceted by region ----
rate_pct_by_region <- df[mean_mort_rate > 0, .(pct = mean(mean_mort_rate <= hotspot_rate_thresh, na.rm = TRUE) * 100),
                         by = region]
rate_pct_by_region[, label := paste0(round(pct), "th percentile")]

g_hist_rate_group_region <- ggplot(df[mean_mort_rate > 0], aes(x = mean_mort_rate, fill = group, color = group)) +
  geom_histogram(alpha = 0.25, bins = 30, position = "identity") +
  geom_vline(xintercept = hotspot_rate_thresh, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_text(data = rate_pct_by_region, inherit.aes = FALSE,
            aes(x = hotspot_rate_thresh, y = Inf, label = label),
            hjust = 1.1, vjust = 1.2, size = 4, angle = 90, color = "grey30") +
  scale_x_sqrt(breaks = c(0, 1, 5, 10, 25, 50)) +
  scale_y_continuous(labels = function(x) x / 100) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_color_manual(values = colorspace::darken(group_colors, 0.4), guide = "none") +
  facet_wrap(~ region, nrow = 1) +
  labs(x = x_label_rate, y = expression("Forest area (km"^2*")")) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(), legend.position = "bottom")
ggsave(file.path(fig_dir, "hist_meanMortRate_overall_byGroup_region.pdf"),
       g_hist_rate_group_region, width = 14, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "hist_meanMortRate_overall_byGroup_region.png"),
       g_hist_rate_group_region, width = 14, height = 4.5, dpi = 300)
plot(g_hist_rate_group_region)

