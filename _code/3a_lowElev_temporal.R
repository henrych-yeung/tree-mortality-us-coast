# LOW ELEVATION TEMPORAL ANALYSIS (0–5 m) — mean mortality rate per year
# Faceted by region; Wetland / Non-wetland / All lines, by 0–5 m elevation bin
# Author: Henry CH Yeung, UVA

library(arrow)
library(data.table)
library(ggplot2)
library(scales)

set.seed(123)

base_dir <- "/Volumes/Henry2/ghostForest_us/_publish_product"   # EDIT to your local path
setwd(base_dir)

CONF_THRESH_STR <- "conf0d5"   # selects which parquet version (and tags outputs)

pq_dir  <- "_full_data/_full"
fig_dir <- paste0("_output/3_temporal_", CONF_THRESH_STR)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

year_cols    <- paste0("annualCount_everg_", 2012:2023)
elev_levels  <- c("0–2 m", "2–5 m")
region_order <- c("Lakes", "Pacific", "Gulf", "Atlantic")

# Florida exclusion bbox (EPSG:5070); set exclude_florida = TRUE to drop FL.
exclude_florida <- FALSE
fl_xmin <- 763003.3315;  fl_xmax <- 1629087.4923
fl_ymin <- 225652.6689;  fl_ymax <- 991763.5669

# ---- Load & filter ----
pq_files <- list.files(pq_dir,
                        pattern = paste0("_", CONF_THRESH_STR, "_.*\\.parquet$"),
                        full.names = TRUE)
cat(paste0("Parquet files (", CONF_THRESH_STR, "): ", length(pq_files), "\n"))

df <- rbindlist(lapply(pq_files, function(f) {
  as.data.table(read_parquet(f, col_select = c(
    "x", "y", "region", "huc6", "dem_3dep", "nlcd", "mean_mort_rate",
    "n_obs_det", "insect_mort", all_of(year_cols)
  )))
}))

df <- df[!is.na(mean_mort_rate) & insect_mort %in% c(0, NA)]
df <- df[dem_3dep >= 0 & dem_3dep < 5]
df[, group := fcase(nlcd == 90,              "Wetland",
                    nlcd %in% c(41, 42, 43), "Non-wetland")]
df <- df[!is.na(group)]
df[, group := factor(group, levels = c("Wetland", "Non-wetland"))]
df[, `:=`(
  elev_bin = factor(fcase(dem_3dep < 2, "0–2 m", dem_3dep < 5, "2–5 m"),
                    levels = elev_levels),
  region   = factor(fcase(tolower(region) == "pacific",  "Pacific",
                          tolower(region) == "atlantic", "Atlantic",
                          tolower(region) == "gulf",     "Gulf",
                          tolower(region) == "lakes",    "Lakes"),
                    levels = region_order)
)]
df <- df[!is.na(region)]
if (exclude_florida)
  df <- df[!(x >= fl_xmin & x <= fl_xmax & y >= fl_ymin & y <= fl_ymax)]
df_base <- df          # all filters except n_obs_det threshold (used for sensitivity ribbon)
df      <- df[n_obs_det > 10 | mean_mort_rate == 0]

# ---- Melt to long & aggregate ----
df_long <- melt(df, id.vars = c("region", "huc6", "elev_bin", "group", "dem_3dep"),
                measure.vars = year_cols, variable.name = "yr", value.name = "ann_mort")
df_long[, year := as.integer(sub("annualCount_everg_", "", yr))][, yr := NULL]
df_long <- df_long[!is.na(ann_mort)]

annual_summary <- df_long[, .(y_mean = mean(ann_mort, na.rm = TRUE),
                               y_sd   = sd(ann_mort,   na.rm = TRUE),
                               n      = .N),
                          by = .(year, region, elev_bin, group)]
annual_summary[, `:=`(se = y_sd / sqrt(n),
                       y_lo = y_mean - 1.96 * y_sd / sqrt(n),
                       y_hi = y_mean + 1.96 * y_sd / sqrt(n))]

# Flag years with < 10% of the max n for that region x elev_bin x group
annual_summary[, n_max := max(n), by = .(region, elev_bin, group)]
annual_summary[, plot_pt := n >= 0.1 * n_max]

y_label <- expression(Mortality ~ rate ~ (ha^{-1} ~ yr^{-1}))

# ================================================================
# Fig: Wetland vs Non-wetland vs All (0–5 m, both elev bins combined)
# ================================================================
ann_group <- df_long[, .(y_mean = mean(ann_mort, na.rm = TRUE),
                          y_sd   = sd(ann_mort,   na.rm = TRUE),
                          n      = .N),
                     by = .(year, region, group)]
ann_group[, `:=`(se = y_sd / sqrt(n), y_lo = y_mean - 1.96 * y_sd / sqrt(n),
                  y_hi = y_mean + 1.96 * y_sd / sqrt(n))]
ann_group[, n_max  := max(n), by = .(region, group)]
ann_group[, plot_pt := n >= 0.1 * n_max]
ann_group[, aggregation := "by_group"]

# "All" line: both groups combined, per year x region
ann_all <- df_long[, .(y_mean = mean(ann_mort, na.rm = TRUE),
                        y_sd   = sd(ann_mort,   na.rm = TRUE),
                        n      = .N),
                   by = .(year, region)]
ann_all[, group := "All"]
ann_all[, `:=`(se = y_sd / sqrt(n), y_lo = y_mean - 1.96 * y_sd / sqrt(n),
               y_hi = y_mean + 1.96 * y_sd / sqrt(n))]
ann_all[, n_max  := max(n), by = .(region)]
ann_all[, plot_pt := n >= 0.1 * n_max]
ann_all[, aggregation := "by_group"]

ann_group <- rbindlist(list(ann_group, ann_all), fill = TRUE)
ann_group[, group := factor(group, levels = c("All", "Wetland", "Non-wetland"))]

colors_group <- c("All" = "grey77", "Wetland" = "#2a9cd5", "Non-wetland" = "#ba8a50", "Hotspot" = "#e96a80")

# ---- Optional background bars: "hotspot" (% pixels > 5 trees/ha/yr), "area"
#      (pixel availability relative to best-coverage year), or "none" ----
bar_mode  <- "none"   # "hotspot" | "area" | "none"
y_mort_hi <- 10        # matches coord_cartesian upper limit

# Restrict bars to (year, region) combos where the line plot has a valid point
valid_yr_reg <- unique(ann_group[plot_pt == TRUE, .(year, region)])

if (bar_mode == "hotspot") {
  bar_df <- df_long[!is.na(region) & region %in% region_order,
                    .(bar_val = mean(ann_mort > 5, na.rm = TRUE)), by = .(year, region)]
  bar_legend <- "Hotspot"; bar_y_label <- "Hotspot proportion (%)"; fname_tag <- "hotspot"
} else if (bar_mode == "area") {
  area_df <- df_long[!is.na(region) & region %in% region_order,
                     .(n_avail = .N), by = .(year, region)]
  area_df[, n_max_region := max(n_avail), by = region]
  area_df[, bar_val := n_avail / n_max_region]
  bar_df <- area_df[, .(year, region, bar_val)]
  bar_legend <- "Area"; bar_y_label <- "Pixel availability (%)"; fname_tag <- "area"
} else if (bar_mode == "none") {
  bar_df <- NULL; fname_tag <- "none"
} else {
  stop("bar_mode must be 'hotspot', 'area', or 'none'")
}

if (!is.null(bar_df)) {
  bar_df <- merge(bar_df, valid_yr_reg, by = c("year", "region"))
  fwrite(bar_df[, .(year, region = tolower(as.character(region)), bar_val)],
         file.path(fig_dir, paste0("bar_", fname_tag, "_by_region_year.csv")))
  a_bar <- y_mort_hi / 1.0   # bar fills the full y range (val is 0-1)
  bar_df[, bar_y := bar_val * a_bar]
}

g_group <- ggplot(ann_group[plot_pt == TRUE], aes(x = year, y = y_mean, color = group, fill = group))

if (!is.null(bar_df)) {
  g_group <- g_group +
    geom_rect(data = bar_df,
              aes(xmin = year - .4, xmax = year + .4, ymin = 0, ymax = bar_y, fill = bar_legend),
              inherit.aes = FALSE, alpha = 0.5)
}

g_group <- g_group +
  # "All" drawn first so it sits underneath the Wetland / Non-wetland lines
  geom_line( data = ann_group[plot_pt == TRUE & group == "All"], linewidth = 0.9) +
  geom_point(data = ann_group[plot_pt == TRUE & group == "All"], size = 2) +
  geom_line( data = ann_group[plot_pt == TRUE & group != "All"], linewidth = 0.9) +
  geom_point(data = ann_group[plot_pt == TRUE & group != "All"], size = 2) +
  scale_color_manual(values = colors_group, name = NULL) +
  scale_x_continuous(breaks = seq(2014, 2023, by = 4))

if (!is.null(bar_df)) {
  g_group <- g_group +
    scale_fill_manual(values = c(colors_group, setNames("#e96a80", bar_legend)),
                      breaks = bar_legend, name = NULL) +
    scale_y_continuous(
      breaks = pretty_breaks(n = 4),
      sec.axis = sec_axis(~ . / a_bar * 100, name = bar_y_label,
                          breaks = pretty(c(0, max(bar_df$bar_val, na.rm = TRUE) * 100), n = 2))
    )
} else {
  g_group <- g_group +
    scale_fill_manual(values = colors_group, guide = "none") +
    scale_y_continuous(breaks = pretty_breaks(n = 4))
}

g_group <- g_group +
  coord_cartesian(ylim = c(0, 10)) +
  facet_wrap(~ region, nrow = 1, axes = "all") +
  labs(x = "Year", y = y_label) +
  theme_classic(base_size = 16) +
  theme(strip.background   = element_blank(),
        panel.grid.major.x = element_line(color = "grey88", linewidth = 0.35),
        strip.text         = element_text(size = 14),
        legend.position    = "top",
        legend.title       = element_blank(),
        legend.background  = element_rect(fill = "transparent", color = NA),
        legend.key         = element_rect(fill = "transparent", color = NA),
        plot.margin        = margin(5, 5, 5, 5))
plot(g_group)

ggsave(file.path(fig_dir, paste0("temporal_mort_wetland_vs_nonwetland_", fname_tag, ".pdf")),
       g_group, width = 15, height = 4.2, dpi = 300)
ggsave(file.path(fig_dir, paste0("temporal_mort_wetland_vs_nonwetland_", fname_tag, ".png")),
       g_group, width = 15, height = 4.2, dpi = 300)

# Annual mortality by region x year x group — feeds 3b's lag-correlation table
fwrite(ann_group, file.path(fig_dir, "stats_temporal.csv"))

# Save the plot object so 3d can stack it above the driver overview figure
saveRDS(g_group, file.path(fig_dir, paste0("g_group_", fname_tag, ".rds")))
