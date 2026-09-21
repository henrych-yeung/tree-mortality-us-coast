# TEMPORAL DRIVER PLOT
# Left y: water level (m); Bars: mutually exclusive stress (drought / hurricane / compound)
# All panels share the same y range; Lakes labels x10 via facetted_pos_scales
# Author: Henry CH Yeung, UVA

library(data.table); library(ggplot2); library(scales); library(ggh4x)

base_dir <- "/Volumes/Henry2/ghostForest_us/_publish_product"   # EDIT to your local path
setwd(base_dir)

CONF_THRESH_STR <- "conf0d5"   # must match the value used in 3b

data_dir <- paste0("_output/3_temporal_", CONF_THRESH_STR)
fig_dir  <- data_dir
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

region_order  <- c("lakes", "pacific", "gulf", "atlantic")
region_labels <- c(lakes = "Lakes", pacific = "Pacific", gulf = "Gulf", atlantic = "Atlantic")

# ---- 1. Load & prepare ----
drv <- fread(file.path(data_dir, "driver_stats_by_region_year.csv"))
cpd <- fread(file.path(data_dir, "compound_stress_by_region_year.csv"))

fct <- function(dt) { dt[, region := factor(region, levels = region_order, labels = region_labels)]; dt }
wl_dt   <- fct(drv[variable == "water_level" & region %in% region_order])
spei_dt <- fct(drv[variable == "spei_anom"   & region %in% region_order])

# Lakes WL / 10 so it shares a comparable axis with coastal gauges;
# Lakes panel labels corrected x10 via facetted_pos_scales
wl_dt[region == "Lakes", c("mean", "ci_lo", "ci_hi") := .(mean / 10, ci_lo / 10, ci_hi / 10)]

# ---- 2. Global normalisation (SPEI + hurr bars onto WL axis) ----
wl_lo  <- min(wl_dt$ci_lo, na.rm = TRUE)
wl_hi  <- max(wl_dt$ci_hi, na.rm = TRUE)
wl_rng <- wl_hi - wl_lo

# ---- 3. Mutually exclusive stress bars (% area) ----
stress_cats   <- c("spei_only", "hurr_only", "compound")
stress_labels <- c(spei_only = "Drought", hurr_only = "Hurricane", compound = "Compound (both)")
stress_colors <- c("Drought" = "#e99b26", "Hurricane" = "#9cadbc", "Compound (both)" = "#e366b3")

pct_cols <- paste0("pct_", stress_cats)
stress_long <- fct(melt(
  cpd[region %in% region_order, c("region", "year", pct_cols), with = FALSE],
  id.vars = c("region", "year"), measure.vars = pct_cols, variable.name = "category", value.name = "pct"
))
stress_long[, category := factor(sub("pct_", "", category), levels = stress_cats, labels = stress_labels)]
stress_long[is.na(pct), pct := 0]

# Scale % onto WL axis (bars use up to half the WL range)
pct_max <- max(stress_long[, .(tot = sum(pct)), by = .(region, year)]$tot, 0.01)
a_pct   <- 0.5 * wl_rng / pct_max

# Stacked ymin/ymax (spei_only at bottom, compound at top)
stress_long <- stress_long[order(region, year, category)]
stress_long[, seg_h := pct * a_pct]
stress_long[, c("ymin", "ymax") := {
  ym <- wl_lo + cumsum(seg_h); list(ym - seg_h, ym)
}, by = .(region, year)]

# ---- 4. Per-panel y scales: same limits everywhere, Lakes labels x10 ----
# wl_lo/wl_hi/a_pct/pct_br are passed as arguments (not looked up from the
# global env) so the resulting formulas stay self-contained after the plot is
# saveRDS()'d and reloaded in a fresh session (see 3d_align_temporal_figures.R).
pct_br <- pretty(c(0, pct_max * 100), n = 2)
make_y <- function(wl_lo, wl_hi, a_pct, pct_br, times10 = FALSE) {
  scale_y_continuous(
    limits = c(wl_lo, wl_hi), name = "Water level (m)",
    labels = if (times10) \(x) x * 10 else waiver(),
    sec.axis = sec_axis(~ (. - wl_lo) / a_pct * 100, name = "Area impacted (%)", breaks = pct_br)
  )
}
facet_y <- list(
  make_y(wl_lo, wl_hi, a_pct, pct_br),
  make_y(wl_lo, wl_hi, a_pct, pct_br),
  make_y(wl_lo, wl_hi, a_pct, pct_br),
  make_y(wl_lo, wl_hi, a_pct, pct_br, times10 = TRUE)
)

# ---- 5. Plot ----
g <- ggplot() +
  geom_rect(data = stress_long[pct > 0],
            aes(xmin = year - .4, xmax = year + .4, ymin = ymin, ymax = ymax, fill = category),
            alpha = 0.85) +
  geom_ribbon(data = wl_dt, aes(x = year, ymin = ci_lo, ymax = ci_hi), fill = "#394955", alpha = 0.2) +
  geom_line(data = wl_dt, aes(x = year, y = mean), color = "#394955", linewidth = 0.9) +
  scale_fill_manual(values = stress_colors, name = "Area impacted by ") +
  scale_x_continuous(breaks = seq(2014, 2023, by = 4)) +
  scale_y_continuous(name = "Water level (m)") +
  coord_cartesian(clip = "off") +
  facet_wrap2(~ region, nrow = 1, axes = "all", scales = "free_y") +
  facetted_pos_scales(y = facet_y) +
  labs(x = "Year") +
  theme_classic(base_size = 16) +
  theme(strip.background   = element_blank(),
        panel.grid.major.x = element_line(color = "grey88", linewidth = 0.35),
        panel.spacing      = unit(0.6, "lines"),
        legend.background  = element_rect(fill = "transparent", color = NA),
        legend.key         = element_rect(fill = "transparent", color = NA),
        legend.position    = "bottom",
        plot.margin        = margin(5, 5, 5, 5))

ggsave(file.path(fig_dir, "temporal_drivers_overview.pdf"), g, width = 18, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "temporal_drivers_overview.png"), g, width = 18, height = 5, dpi = 300)
plot(g)

# Save the plot object so 3d can stack it below the group-comparison figure (3a)
saveRDS(g, file.path(fig_dir, "g_driver.rds"))
cat("Saved to:", fig_dir, "\n")
