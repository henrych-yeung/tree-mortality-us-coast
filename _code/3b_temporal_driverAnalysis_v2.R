# TEMPORAL DRIVER ANALYSIS v2 — mutually exclusive stress categories
# Predictors : pct_spei_only + pct_hurr_only + pct_compound + water_level
# Lakes      : pct_spei_only + water_level (no hurricane data)
# Subgroups  : All / Wetland / Non-wetland, 0-5 m elevation
# Analysis   : Pearson correlation with lags 0-2 (bivariate, exploratory)
#
# MLR was considered but dropped: (1) n <= 12 years/region leaves too few
# degrees of freedom for 4 predictors; (2) the mutually exclusive stress
# fractions are collinear by construction; (3) mortality rate and hotspot
# fraction both violate OLS assumptions. A mixed-effects model on HUC6-level
# data with region as a random effect would give valid inference instead.
# Author: Henry CH Yeung, UVA

library(data.table)

base_dir <- "/Volumes/Henry2/ghostForest_us/_publish_product"   # EDIT to your local path
setwd(base_dir)

CONF_THRESH_STR <- "conf0d5"   # must match the value used in scripts 3 & 3a

data_dir <- paste0("_output/3_temporal_", CONF_THRESH_STR)
fig_dir  <- paste0("_output/3_driverAnalysis_", CONF_THRESH_STR)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

region_order <- c("lakes", "pacific", "gulf", "atlantic")

# ================================================================
# 1. Load & prepare drivers
# ================================================================
drv <- fread(file.path(data_dir, "driver_stats_by_region_year.csv"))
cpd <- fread(file.path(data_dir, "compound_stress_by_region_year.csv"))

wl_wide  <- dcast(drv[variable == "water_level"], region + year ~ variable, value.var = "mean")
cpd_wide <- cpd[, .(region, year, pct_spei_only, pct_hurr_only, pct_compound)]

drv_wide <- merge(wl_wide, cpd_wide, by = c("region", "year"), all = TRUE)
# Non-hurricane years / regions -> 0 for hurricane-dependent columns
drv_wide[is.na(pct_hurr_only), pct_hurr_only := 0]
drv_wide[is.na(pct_compound),  pct_compound  := 0]

# ================================================================
# 2. Mortality: All / Wetland / Non-wetland, 0-5 m
# ================================================================
mort_all <- fread(file.path(data_dir, "stats_temporal.csv"))

get_mort <- function(grp) {
  mort_all[aggregation == "by_group" & group == grp & plot_pt == TRUE,
           .(y = weighted.mean(y_mean, n, na.rm = TRUE)),
           by = .(region = tolower(as.character(region)), year)]
}

mort <- list(all = get_mort("All"), wetland = get_mort("Wetland"), nonwetland = get_mort("Non-wetland"))
dat  <- lapply(mort, function(m) merge(m, drv_wide, by = c("region", "year"), all.x = TRUE))

# ================================================================
# 3. Variable selection  <- EDIT HERE
# ================================================================
# Predictors per region. Comment out a line to exclude that predictor.
# In bivariate correlation, collinear predictors don't bias individual r
# estimates; dropping one is a scientific redundancy choice, not a
# statistical requirement.
region_vars <- list(
  pacific  = c("pct_spei_only", "water_level"),   # no Pacific hurricane track
  atlantic = c("pct_spei_only", "pct_hurr_only", "pct_compound", "water_level"),
  gulf     = c("pct_spei_only", "pct_hurr_only", "pct_compound", "water_level"),
  lakes    = c("pct_spei_only", "water_level")    # no Great Lakes hurricane data
)

all_pred_labels <- c(
  pct_spei_only = "Drought only",
  pct_hurr_only = "Hurricane only",
  pct_compound  = "Compound (both)",
  water_level   = "Water level (m)"
)

use_vars   <- unique(unlist(region_vars))   # union across regions
use_labels <- all_pred_labels[use_vars]

scenarios <- list(
  s1 = list(region_vars = region_vars, all_vars = use_vars, labels = use_labels)
)

# ================================================================
# 4. Lag correlation (Pearson, lags 0-2)
# ================================================================
sig_stars <- function(p) fcase(p < 0.001, "***", p < 0.01, "**", p < 0.05, "*", default = "")

run_lag_cors <- function(d_reg, sg_label, sc) {
  rbindlist(lapply(region_order, function(reg) {
    vars <- sc$region_vars[[reg]]
    d    <- d_reg[[reg]][order(year)]
    rbindlist(lapply(vars, function(v) {
      rbindlist(lapply(0:2, function(lag) {
        x_lag <- data.table::shift(d[[v]], lag, type = "lag")
        ok    <- !is.na(d$y) & !is.na(x_lag)
        if (sum(ok) < 4) return(NULL)
        ct <- cor.test(d$y[ok], x_lag[ok], method = "pearson")
        data.table(subgroup = sg_label, region = reg,
                   predictor = factor(v, levels = sc$all_vars),
                   lag = lag, r = round(ct$estimate, 3),
                   p = round(ct$p.value, 4), n = sum(ok),
                   sig = sig_stars(ct$p.value))
      }))
    }))
  }))
}

# ================================================================
# 5. Run for each subgroup x scenario
# ================================================================
subgroup_meta <- list(all = dat$all, wetland = dat$wetland, nonwetland = dat$nonwetland)

all_lag <- list()
for (sg_name in names(subgroup_meta)) {
  d_reg <- lapply(setNames(region_order, region_order), function(reg) subgroup_meta[[sg_name]][region == reg])

  for (sc_name in names(scenarios)) {
    sc  <- scenarios[[sc_name]]
    key <- paste0(sg_name, "_", sc_name)
    all_lag[[key]] <- run_lag_cors(d_reg, key, sc)
  }
}

# ================================================================
# 6. Publishable table — Pearson correlations (all lags), wide CSV
# ================================================================
sg_display <- c(all_s1 = "All", wetland_s1 = "Wetland", nonwetland_s1 = "Non-wetland")
sg_key     <- c("All" = "al", "Wetland" = "wl", "Non-wetland" = "nw")

pub_long <- rbindlist(all_lag)
pub_long[, `:=`(
  Subgroup  = sg_display[subgroup],
  Region    = region,
  Predictor = use_labels[as.character(predictor)],
  cell      = sprintf("%.2f%s", r, sig),
  col_key   = paste0(sg_key[sg_display[subgroup]], "_lag", lag)
)]

pub_wide <- dcast(pub_long, Region + Predictor ~ col_key, value.var = "cell")

# Enforce column order: All lag 0-2, Wetland lag 0-2, Non-wetland lag 0-2
col_order <- c("Region", "Predictor", paste0("al_lag", 0:2), paste0("wl_lag", 0:2), paste0("nw_lag", 0:2))
for (col in col_order[!col_order %in% names(pub_wide)]) pub_wide[, (col) := NA_character_]
setcolorder(pub_wide, col_order)

table_region_order <- c("lakes", "pacific", "gulf", "atlantic")
table_pred_order   <- c("Water level (m)", "Hurricane only", "Drought only", "Compound (both)")
pub_wide[, `:=`(Region = factor(Region, levels = table_region_order),
                Predictor = factor(Predictor, levels = table_pred_order))]
setorder(pub_wide, Region, Predictor)

fwrite(pub_wide, file.path(fig_dir, "lagcor_v2_table.csv"))
cat("Saved to:", fig_dir, "\n")
