# TEMPORAL DRIVER PREP v2 — pixel-level stats from filtered forest pixels
# Pixel filter: same as 1_spatialAnalysis.R (0–5 m elev, NLCD 41/42/43/90, insect filter)
# SPEI & hurricane: raster extraction at filtered pixel locations
# Water level: NOAA gauges
# Author: Henry CH Yeung, UVA

library(data.table); library(terra); library(sf); library(arrow)
library(ggplot2); library(scales); library(patchwork)

base_dir <- "/Volumes/Henry2/ghostForest_us/_publish_product"   # EDIT to your local path
setwd(base_dir)

CONF_THRESH_STR <- "conf0d5"   # selects parquet version (and tags outputs)

pq_dir       <- "_driver_product/_full"
aux_dir      <- "_driver_product/_aux_temporal"    # spei / hurricane / water level / region
out_dir      <- paste0("_output/3_temporal_", CONF_THRESH_STR)
driver_years <- 2012:2023
region_order  <- c("lakes", "pacific", "gulf", "atlantic")
region_labels <- c(lakes = "Lakes", pacific = "Pacific", gulf = "Gulf", atlantic = "Atlantic")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ================================================================
# 1. Load & filter forest pixels from parquet (filtered by conf)
# ================================================================
pq_files <- list.files(pq_dir, pattern = paste0("_", CONF_THRESH_STR, "_.*\\.parquet$"), full.names = TRUE)
cat(paste0("Parquet files (", CONF_THRESH_STR, "): ", length(pq_files), "\n"))

df <- rbindlist(lapply(pq_files, function(f) {
  as.data.table(read_parquet(f, col_select = c(
    "x", "y", "region", "huc6", "dem_3dep", "nlcd",
    "mean_mort_rate", "n_obs_det", "insect_mort"
  )))
}))

df <- df[!is.na(mean_mort_rate) & insect_mort %in% c(0, NA)]
df <- df[dem_3dep >= 0 & dem_3dep < 5]
df[, group := fcase(nlcd == 90L,               "Wetland",
                    nlcd %in% c(41L, 42L, 43L), "Non-wetland")]
df <- df[!is.na(group)]
df[, region := tolower(region)]
df <- df[region %in% region_order]
df <- df[n_obs_det > 10 | mean_mort_rate == 0]

# Deduplicate (a pixel can appear in multiple parquet files)
coords_dt <- unique(df, by = c("x", "y"))
cat("Forest pixels for driver sampling:", nrow(coords_dt), "\n")

n_px_reg <- coords_dt[, .(n_pixels = .N), by = region]   # fixed denominator per region

pts_alb <- vect(coords_dt[, .(x, y)], geom = c("x", "y"), crs = "EPSG:5070")   # parquet's native CRS

# ================================================================
# 2. Hurricane EF-scale damage raster — extract at pixel locations
# ================================================================
dmg_rast  <- rast(file.path(aux_dir, "accum_damage_simpleSum_byYear_2012_2023.tif"))
pts_wgs   <- project(pts_alb, crs(dmg_rast))
dmg_years <- as.integer(sub("acc_", "", names(dmg_rast)))   # 2012, 2014-2023 (no 2013)

dmg_ex <- as.data.table(terra::extract(dmg_rast, pts_wgs))[, -"ID"]
setnames(dmg_ex, as.character(dmg_years))
dmg_ex[is.na(dmg_ex)] <- 0   # outside model domain = no damage
dmg_ex[, pixel_id := .I]

dmg_long <- melt(dmg_ex, id.vars = "pixel_id", measure.vars = as.character(dmg_years),
                 variable.name = "year", value.name = "ef_damage")
dmg_long[, `:=`(year = as.integer(as.character(year)), region = coords_dt$region[pixel_id])]

# ================================================================
# 3. SPEI raster — extract at pixel locations (already Albers)
# ================================================================
spei_rast <- rast(file.path(aux_dir, "spei1y_yearly_2012_2023.tif"))
names(spei_rast) <- driver_years

spei_ex <- as.data.table(terra::extract(spei_rast, pts_alb))[, -"ID"]
setnames(spei_ex, as.character(driver_years))
spei_ex[, pixel_id := .I]

spei_long <- melt(spei_ex, id.vars = "pixel_id", measure.vars = as.character(driver_years),
                  variable.name = "year", value.name = "spei_anom")
spei_long[, `:=`(year = as.integer(as.character(year)), region = coords_dt$region[pixel_id])]

# ================================================================
# 4. Compound stress (mutually exclusive) — pixel-level; denominator =
#    filtered forest pixels. Categories: drought only, hurricane only,
#    compound (both); these are the only drivers plotted downstream.
# ================================================================
spei_neg1_region <- spei_long[!is.na(spei_anom), .(n_spei_neg1 = sum(spei_anom < -1)), by = .(region, year)]

# Hurricane-dependent counts (inner join -> sparse for Lakes / 2013)
dmg_spei_px <- merge(dmg_long[, .(pixel_id, year, region, ef_damage)],
                     spei_long[!is.na(spei_anom), .(pixel_id, year, spei_anom)],
                     by = c("pixel_id", "year"))

compound_hurr <- dmg_spei_px[, .(
  n_spei_only = sum(spei_anom < -1 & ef_damage == 0),
  n_hurr_only = sum(ef_damage >= 1 & spei_anom > -1),
  n_compound  = sum(ef_damage >= 1 & spei_anom < -1)
), by = .(region, year)]

# Full region x year grid so Lakes / 2013 get 0 for hurricane columns
full_grid <- CJ(region = region_order, year = driver_years)
compound_region <- merge(full_grid, compound_hurr, by = c("region", "year"), all.x = TRUE)
hurr_cols <- c("n_hurr_only", "n_compound")
for (col in hurr_cols) compound_region[is.na(get(col)), (col) := 0L]

# Attach SPEI-only counts (covers all regions x years)
compound_region <- merge(compound_region, spei_neg1_region, by = c("region", "year"), all.x = TRUE)
# For non-hurricane rows n_spei_only is NA -> all SPEI<-1 pixels are "spei only"
compound_region[is.na(n_spei_only), n_spei_only := n_spei_neg1]
compound_region[, n_spei_neg1 := NULL]

compound_region <- merge(compound_region, n_px_reg, by = "region")
pct_vars <- c("n_spei_only", "n_hurr_only", "n_compound")
for (v in pct_vars)
  compound_region[, paste0("pct_", sub("n_", "", v)) := get(v) / n_pixels]

# ================================================================
# 5. Water level — NOAA gauges
# ================================================================
regions_sf <- st_read(file.path(aux_dir, "region.gpkg"), quiet = TRUE)
regions_sf$region <- tolower(regions_sf$region)

wl    <- st_read(file.path(aux_dir, "waterLevel_noaa_msl_summary.gpkg"), quiet = TRUE)
wl_dt <- as.data.table(st_drop_geometry(
  st_join(st_transform(wl, st_crs(regions_sf)), regions_sf, join = st_within)
))[!is.na(region)]

mean_cols <- paste0("mean_", driver_years)
wl_long   <- melt(wl_dt[, c("station_id", "region", ..mean_cols)],
                  id.vars = c("station_id", "region"), measure.vars = mean_cols,
                  variable.name = "yr", value.name = "wl_val")
wl_long[, year := as.integer(sub("mean_", "", yr))][, yr := NULL]

stats_wl <- wl_long[!is.na(wl_val), .(
  variable = "water_level", mean = mean(wl_val), se = sd(wl_val) / sqrt(.N)
), by = .(region, year)]

# ================================================================
# 6. Export CSVs
# ================================================================
driver_stats <- copy(stats_wl)
driver_stats[, `:=`(ci_lo = mean - 1.96 * se, ci_hi = mean + 1.96 * se)]
fwrite(driver_stats, file.path(out_dir, "driver_stats_by_region_year.csv"))

fwrite(compound_region[region %in% region_order],
       file.path(out_dir, "compound_stress_by_region_year.csv"))
cat("CSVs saved to:", out_dir, "\n")
