# LANDCOVER GROUP ANALYSIS (Non-wetland vs Wetland) — bootstrap aggregation
# Runs propensity-score matching + line plots on all 25 bootstrap samples.
# Bar plot CI = bootstrap mean +/- 1.96 x SD across bootstrap estimates.
# Line plots = loess/lm fit on pooled bin-level summaries (not raw data).
# Author: Henry CH Yeung, UVA

library(data.table)
library(ggplot2)
library(tidyverse)
library(scales)
library(MatchIt)
library(marginaleffects)
library(colorspace)
library(corrplot)
library(ggpubr)

# Resolve potential conflicts: some transitive deps (Hmisc/plyr) mask dplyr::summarize
summarize <- dplyr::summarize
summarise <- dplyr::summarise
filter    <- dplyr::filter
select    <- dplyr::select

seed <- 123
set.seed(seed)

base_dir <- "/Volumes/Henry2/ghostForest_us/_publish_product"   # EDIT to your local path
setwd(base_dir)

CONF_THRESH_STR <- "conf0d5"   # selects which bootstrap version to read

bs_dir <- paste0("_driver_product/_bs_highElev_", CONF_THRESH_STR, "_0d02_20000")
bs_files <- sort(list.files(bs_dir, pattern = "sample_bs_highElev_b\\d+\\.csv$", full.names = TRUE))
cat(paste0("Bootstrap files found: ", length(bs_files), " in ", bs_dir, "\n"))

# Output directories (tagged with conf threshold for traceability)
output_directory     <- paste0("_output/4_landcoverAnalysis_", CONF_THRESH_STR, "/data")
output_directory_fig <- paste0("_output/4_landcoverAnalysis_", CONF_THRESH_STR, "/fig")
if (!dir.exists(output_directory))     dir.create(output_directory,     recursive = TRUE)
if (!dir.exists(output_directory_fig)) dir.create(output_directory_fig, recursive = TRUE)

# Region subset: NULL = all regions, or a single value ("pacific"/"gulf"/
# "atlantic"/"lakes") to filter. When a single region is set, "region" is
# dropped from `controls` to avoid a 1-level-factor contrast error.
region_filter <- c()

# Florida exclusion bbox (EPSG:5070); set exclude_florida = TRUE to drop FL.
exclude_florida <- FALSE
fl_xmin <- 763003.3315;  fl_xmax <- 1629087.4923
fl_ymin <- 225652.6689;  fl_ymax <- 991763.5669

# Controls (used in matching formula and outcome model)
controls <- c(
  "tree_cover", "forest_type", "ppt",
  "dist_flowline", "density_flowline",
  "wet_freq_z", "wet_sev_z", "wet_freq", "tmax_z", "spei_freq",
  "insect_defo", "insect_mort",
  "hurr_ws", "region", "tpa_live"
)
if (!is.null(region_filter)) controls <- setdiff(controls, "region")

palette_lc   <- c("Non-wetland" = "#ba8a50", "Wetland" = "#2a9cd5")
y_label      <- expression(Mortality ~ rate ~ (ha^-1 ~ yr^-1))
scatter_vars <- c("ppt","wet_sev_z")
no_clip_vars <- c()  # discrete: skip percentile clipping in get_bins

# ---- Load & prepare one bootstrap file ----
prepare_data <- function(bs_file) {
  df <- fread(bs_file)
  df <- df[nlcd %in% c(41, 42, 43, 90)]
  df <- df[hurr_ws == 0]   # ==0 non-storm-impacted; >0 storm-impacted zone only
  df[, lc_group := factor(fifelse(nlcd %in% c(41, 42, 43), "Non-wetland", "Wetland"),
                          levels = c("Wetland", "Non-wetland"))]
  df <- df[!is.na(mean_mort_rate) & (n_obs_det > 10 | mean_mort_rate == 0)]
  if (!is.null(region_filter)) df <- df[region == region_filter]
  if (exclude_florida)
    df <- df[!(x >= fl_xmin & x <= fl_xmax & y >= fl_ymin & y <= fl_ymax)]
  pa_cols <- c("insect_defo", "insect_mort", "gmw", "sss", "tidal_range")
  df[, (pa_cols) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = pa_cols]
  df[, r95p := r95p * 365.25]
  df[, fd := fd * 365.25]
  df[, tn10p := tn10p * 365.25]
  df
}

#### PSM — returns marginal predictions only ####
run_psm <- function(data, case_id, controls) {

  cat(paste0("\n=== ", case_id, " ===\n"))

  keep_cols <- c("mean_mort_rate", "lc_group", controls)
  pre_data  <- data[, ..keep_cols] %>% na.omit() %>% as.data.frame()
  pre_data$treat <- as.integer(pre_data$lc_group == "Wetland")
  cat(paste0(nrow(pre_data), " rows after NA removal\n"))

  m.out <- matchit(as.formula(paste("treat ~", paste(controls, collapse = " + "))),
                   data     = pre_data,
                   estimand = "ATT",
                   replace  = FALSE,
                   method   = "nearest",
                   distance = "glm",
                   caliper  = 0.01,
                   verbose  = FALSE)

  m.data <- match.data(m.out)
  m_summary <- summary(m.out)

  # Standardized mean differences for the love plot (pre-match "All" vs "Matched").
  # summary.matchit stores these in sum.all / sum.matched under "Std. Mean Diff.".
  smd_col <- "Std. Mean Diff."
  smd_df  <- data.frame(
    covariate   = rownames(m_summary$sum.all),
    smd_all     = abs(m_summary$sum.all[,     smd_col]),
    smd_matched = abs(m_summary$sum.matched[, smd_col]),
    bs_id       = case_id,
    row.names   = NULL
  )

  fit <- glm(as.formula(paste("mean_mort_rate ~ treat *(",
                              paste(controls, collapse = " + "), ")")),
             family  = Gamma(link = "log"),
             data    = m.data,
             weights = weights,
             control = glm.control(maxit = 100))

  avg_pred <- avg_predictions(fit,
                              variables  = "treat",
                              vcov       = ~subclass,
                              newdata    = subset(m.data, treat == 1),
                              wts        = "weights",
                              conf_level = 0.95)
  marg_df          <- as.data.frame(avg_pred)
  marg_df$treat    <- as.integer(as.character(marg_df$treat))
  marg_df$lc_group <- factor(ifelse(marg_df$treat == 1, "Wetland", "Non-wetland"),
                             levels = c("Wetland", "Non-wetland"))
  marg_df$bs_id    <- case_id

  # GLM fit statistics
  fit_s    <- summary(fit)
  fit_stat <- data.frame(
    bs_id              = case_id,
    n_matched          = nrow(m.data),
    null_deviance      = fit_s$null.deviance,
    residual_deviance  = fit_s$deviance,
    deviance_explained = round(1 - fit_s$deviance / fit_s$null.deviance, 4),
    aic                = fit_s$aic,
    dispersion         = fit_s$dispersion
  )

  m.data$lc_group <- factor(ifelse(m.data$treat == 1, "Wetland", "Non-wetland"),
                            levels = c("Wetland", "Non-wetland"))

  # Pixel counts: per group and total (pre-match and matched)
  n_pre     <- pre_data %>% count(lc_group) %>% rename(n_pre = n)
  n_matched <- m.data   %>% count(lc_group) %>% rename(n_matched = n)
  n_counts  <- left_join(n_pre, n_matched, by = "lc_group") %>%
    mutate(n_total_pre     = sum(n_pre),
           n_total_matched = sum(n_matched),
           bs_id           = case_id)

  marg_df <- left_join(marg_df,
                       n_counts %>% select(lc_group, n_matched, n_total_matched),
                       by = "lc_group")

  raw_means <- m.data %>%
    group_by(lc_group) %>%
    summarize(raw_mean = mean(mean_mort_rate, na.rm = TRUE), .groups = "drop") %>%
    left_join(n_counts %>% select(lc_group, n_pre, n_matched, n_total_pre, n_total_matched, bs_id),
              by = "lc_group") %>%
    mutate(bs_id = case_id)

  list(marginal = marg_df, matched = m.data, fit_stat = fit_stat,
       raw_means = raw_means, smd = smd_df)
}

cut_lo <- 0.025
cut_hi <- 0.975

#### Bin-level summaries for line plots ####
get_bins <- function(matched_data, vars, case_id, n_bins = 15, no_clip_vars = NULL,
                     y_var = "mean_mort_rate") {
  all_bins <- list()

  for (var in vars) {

    # Crop to 2.5th-97.5th percentile (skip for count/binary vars in no_clip_vars)
    cropped <- if (!is.null(no_clip_vars) && var %in% no_clip_vars) {
      matched_data
    } else {
      p.lower <- quantile(matched_data[[var]], cut_lo, na.rm = TRUE)
      p.upper <- quantile(matched_data[[var]], cut_hi, na.rm = TRUE)
      matched_data %>% dplyr::filter(.data[[var]] >= p.lower, .data[[var]] <= p.upper)
    }

    bin_df <- cropped %>%
      mutate(bin = cut(.data[[var]], breaks = n_bins, labels = FALSE)) %>%
      group_by(lc_group, bin) %>%
      summarize(
        x_mean = mean(.data[[var]],    na.rm = TRUE),
        y_mean = mean(.data[[y_var]],  na.rm = TRUE),
        y_sd   = sd(.data[[y_var]],    na.rm = TRUE),
        y_se   = sd(.data[[y_var]],    na.rm = TRUE) / sqrt(n()),
        n      = n(),
        .groups = "drop"
      ) %>%
      mutate(variable = var, bs_id = case_id)

    all_bins[[var]] <- bin_df
  }

  do.call(rbind, all_bins)
}

#### Main loop over bootstrap files ####
all_marginals <- list()
all_bins      <- list()
all_fit_stats <- list()
all_raw_means <- list()
all_matched   <- list()
all_smd       <- list()

for (bs_file in bs_files) {
  bs_id <- sub("sample_bs_highElev_", "", tools::file_path_sans_ext(basename(bs_file)))
  df    <- prepare_data(bs_file)

  res                    <- run_psm(df, case_id = bs_id, controls = controls)
  all_marginals[[bs_id]] <- res$marginal
  all_fit_stats[[bs_id]] <- res$fit_stat
  all_raw_means[[bs_id]] <- res$raw_means
  all_matched[[bs_id]]   <- res$matched
  all_smd[[bs_id]]       <- res$smd

  all_bins[[bs_id]] <- get_bins(res$matched, vars = scatter_vars, case_id = bs_id,
                                n_bins = 15, no_clip_vars = no_clip_vars)
}

marginals_all <- do.call(rbind, all_marginals)
bins_all      <- as.data.frame(dplyr::bind_rows(all_bins))
fit_stats_all <- do.call(rbind, all_fit_stats)
raw_means_all <- do.call(rbind, all_raw_means)
smd_all_bs    <- do.call(rbind, all_smd)

fwrite(marginals_all, file.path(output_directory, "bs_marginal_all.csv"))
fwrite(bins_all,      file.path(output_directory, "bs_bins_all.csv"))
fwrite(fit_stats_all, file.path(output_directory, "bs_glm_fit_stats.csv"))
fwrite(smd_all_bs,    file.path(output_directory, "bs_balance_smd_all.csv"))
print(fit_stats_all)

#### Covariate balance love plot: bootstrap mean +/- SD of |SMD| ####
# Mirrors plot(summary(m.out)): open circle = pre-match ("All"), filled =
# "Matched". Points = mean across the 25 bootstraps; horizontal bars = +/- 1
# SD. Reference lines at 0 (solid), 0.05 (dashed), 0.1 (solid).
agg_smd <- smd_all_bs %>%
  group_by(covariate) %>%
  summarize(
    all_mean     = mean(smd_all,     na.rm = TRUE),
    all_sd       = sd(smd_all,       na.rm = TRUE),
    matched_mean = mean(smd_matched, na.rm = TRUE),
    matched_sd   = sd(smd_matched,   na.rm = TRUE),
    .groups = "drop"
  )
fwrite(agg_smd, file.path(output_directory, "agg_balance_smd.csv"))

# Long format: one row per covariate x sample, with mean and sd columns
smd_long <- agg_smd %>%
  pivot_longer(
    cols          = c(all_mean, all_sd, matched_mean, matched_sd),
    names_to      = c("sample", ".value"),
    names_pattern = "(all|matched)_(mean|sd)"
  ) %>%
  mutate(sample = factor(ifelse(sample == "all", "All", "Matched"), levels = c("All", "Matched")))

# Order covariates so the largest pre-match imbalance sits at the top
cov_order <- agg_smd %>% arrange(all_mean) %>% pull(covariate)
smd_long  <- smd_long %>% mutate(covariate = factor(covariate, levels = cov_order))

g_love <- ggplot(smd_long, aes(x = mean, y = covariate, shape = sample)) +
  geom_vline(xintercept = 0,    color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0.05, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.1,  color = "black", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = pmax(mean - sd, 0), xmax = mean + sd),
                 height = 0, linewidth = 0.6, color = "black") +
  geom_point(size = 2.6, fill = "white", color = "black") +
  scale_shape_manual(values = c("All" = 21, "Matched" = 16)) +
  labs(x = "Absolute Standardized\nMean Difference", y = NULL, shape = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position    = c(0.85, 0.15),
        legend.background  = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"))
plot(g_love)
ggsave(file.path(output_directory_fig, "balance_loveplot_bs_agg.pdf"), g_love, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_directory_fig, "balance_loveplot_bs_agg.png"), g_love, width = 6, height = 5, dpi = 300)

# Correlation among scatter_vars (first bootstrap sample)
matched_first <- all_matched[[1]]
df_cor <- matched_first[, scatter_vars] %>%
  mutate(across(everything(), as.numeric)) %>%
  na.omit()
cor_mat <- cor(df_cor, method = "pearson", use = "pairwise.complete.obs")
print(round(cor_mat, 3))
pdf(file.path(output_directory_fig, "scatter_vars_correlation.pdf"), width = 5, height = 5)
jpeg(file.path(output_directory_fig, "scatter_vars_correlation.jpg"), width = 5, height = 5, units = "in", res = 300)
corrplot(cor_mat, method = "circle", type = "upper",
         addCoef.col = "black", number.cex = 1.2, tl.col = "black", tl.srt = 45, diag = FALSE)
dev.off()
corrplot(cor_mat, method = "circle", type = "upper",
         addCoef.col = "black", number.cex = 1.2, tl.col = "black", tl.srt = 45, diag = FALSE)

#### Helper: bar plot with per-bar value labels + a bracket showing the % diff ####
make_bar_diff_plot <- function(agg_df, percent_diff, y_label_plot = y_label,
                               fill_values = palette_lc) {
  wet_top    <- agg_df$conf.high[agg_df$lc_group == "Wetland"]
  nonwet_top <- agg_df$conf.high[agg_df$lc_group == "Non-wetland"]
  y_max      <- max(agg_df$conf.high)
  label_gap  <- y_max * 0.09
  bracket_y  <- y_max * 1.28

  ggplot(agg_df, aes(x = lc_group, y = estimate, fill = lc_group)) +
    geom_bar(stat = "identity", position = "dodge", color = "black",
             linewidth = 0.8, width = 0.5, alpha = 0.6) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = 0.2, linewidth = 0.8, color = "black") +
    geom_text(aes(y = conf.high + label_gap, label = sprintf("%.1f", estimate)),
              size = 5, color = "black") +
    # Bracket connecting the two bar tops, with the % difference on it
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
             linewidth = 0.6, color = "black") +
    annotate("segment", x = 1, xend = 1, y = wet_top + label_gap * 2.2, yend = bracket_y,
             linewidth = 0.6, color = "black") +
    annotate("segment", x = 2, xend = 2, y = nonwet_top + label_gap * 2.2, yend = bracket_y,
             linewidth = 0.6, color = "black") +
    annotate("text", x = 1.5, y = bracket_y + label_gap * 1.5,
             label = sprintf("%+.1f%%", percent_diff), size = 6, color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = scales::breaks_extended(n = 4)) +
    scale_fill_manual(values = fill_values) +
    labs(x = NULL, y = y_label_plot) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none")
}

#### Raw matched mean: bootstrap mean +/- 1.96 * SD across bootstrap samples ####
agg_raw <- raw_means_all %>%
  group_by(lc_group) %>%
  summarize(
    se_bs     = sd(raw_mean),
    estimate  = mean(raw_mean),
    conf.low  = pmax(estimate - 1.96 * se_bs, 0),
    conf.high = estimate + 1.96 * se_bs,
    .groups   = "drop"
  ) %>%
  select(lc_group, conf.low, conf.high, estimate)

percent_diff_raw <- round(
  (agg_raw$estimate[agg_raw$lc_group == "Wetland"] -
     agg_raw$estimate[agg_raw$lc_group == "Non-wetland"]) /
    agg_raw$estimate[agg_raw$lc_group == "Non-wetland"] * 100, 1)
cat(paste0("Raw matched mean — Wetland vs Non-wetland: ", percent_diff_raw, "%\n"))

fwrite(agg_raw, file.path(output_directory, "agg_raw_matched_mean.csv"))

g_bar_raw <- make_bar_diff_plot(agg_raw, percent_diff_raw)
plot(g_bar_raw)
ggsave(file.path(output_directory_fig, "landcover_rawmean_bs_agg.pdf"), g_bar_raw, width = 3.5, height = 4, dpi = 300)
ggsave(file.path(output_directory_fig, "landcover_rawmean_bs_agg.png"), g_bar_raw, width = 3.5, height = 4, dpi = 300)

#### Line plots on pooled bin-level summaries ####
plot_lines_pooled <- function(bins_all, vars, n_bins = 15,
                              point_size_range = c(1, 5),
                              point_size_breaks = NULL,
                              min_n = NULL,
                              x_labels = NULL,
                              y_label_plot = y_label,
                              y_max = NULL,
                              file_tag = "",
                              smooth_method = list(),
                              y_max_by_var  = list(),
                              discrete_vars = NULL) {

  for (var in vars) {
    bin_df <- bins_all[bins_all$variable == var, ]

    # Pool across bootstraps. Continuous vars: group by the per-bootstrap
    # bin index from get_bins. Discrete vars (e.g. insect_mort): group by
    # the rounded x value since bin indices may not align across bootstraps
    # when the data range varies.
    if (!is.null(discrete_vars) && var %in% discrete_vars) {
      bin_df <- bin_df %>% mutate(bin = round(x_mean, 4))
    }
    plot_df <- bin_df %>%
      group_by(lc_group, bin) %>%
      summarize(
        x_mean = mean(x_mean, na.rm = TRUE),
        y_mean = mean(y_mean, na.rm = TRUE),
        y_sd   = sqrt(mean(y_sd^2, na.rm = TRUE)),  # RMS of per-bootstrap within-bin SDs
        y_se   = sqrt(mean(y_se^2, na.rm = TRUE)),
        n      = sum(n),
        .groups = "drop"
      ) %>%
      filter(!is.na(bin))
    if (!is.null(min_n)) plot_df <- plot_df %>% filter(n >= min_n)

    scatter_out <- plot_df %>%
      mutate(variable = var) %>%
      select(variable, lc_group, bin, x_mean, y_mean, y_sd, y_se, n) %>%
      arrange(lc_group, x_mean)
    fwrite(scatter_out,
           file.path(output_directory_fig, paste0("scatter_bin_data_", var, file_tag, ".csv")),
           row.names = FALSE)

    # Fit lm on re-binned plot_df, weighted by 1/y_se^2 — consistent with plotted line
    rebin_stats <- do.call(rbind, lapply(c("Wetland", "Non-wetland"), function(grp) {
      sub_df    <- plot_df[plot_df$lc_group == grp, ]
      fit       <- lm(y_mean ~ x_mean, data = sub_df, weights = 1 / y_se^2)
      coef_s    <- summary(fit)$coefficients
      has_slope <- nrow(coef_s) >= 2
      data.frame(
        variable     = var, lc_group = grp,
        n_bins       = nrow(sub_df),
        n_total      = sum(sub_df$n),
        intercept    = round(coef_s[1, "Estimate"],  5),
        intercept_se = round(coef_s[1, "Std. Error"], 5),
        slope        = if (has_slope) round(coef_s[2, "Estimate"],   5) else NA_real_,
        slope_se     = if (has_slope) round(coef_s[2, "Std. Error"], 5) else NA_real_,
        p_value      = if (has_slope) round(coef_s[2, "Pr(>|t|)"],   5) else NA_real_,
        r2           = round(summary(fit)$r.squared, 4),
        x_min        = min(sub_df$x_mean, na.rm = TRUE),
        x_max        = max(sub_df$x_mean, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    print(rebin_stats)
    fwrite(rebin_stats, file.path(output_directory_fig, paste0("rebin_lm_stats_", var, ".csv")), row.names = FALSE)

    sm <- if (var %in% names(smooth_method)) smooth_method[[var]] else "both"

    # y-axis limit priority: y_max_by_var[[var]] > global y_max > auto
    y_lim <- if (var %in% names(y_max_by_var)) {
      y_max_by_var[[var]]
    } else if (!is.null(y_max)) {
      y_max
    } else {
      max(plot_df$y_mean + plot_df$y_se, na.rm = TRUE) * 1.1
    }

    g_line <- ggplot(plot_df, aes(x = x_mean, y = y_mean, color = lc_group))

    if (sm %in% c("both", "loess"))
      g_line <- g_line +
      geom_smooth(method = "loess", span = 0.6, se = TRUE, linewidth = 1, alpha = 0.4,
                  aes(weight = 1 / y_se^2, fill = after_scale(lighten(color, 0.5))))

    if (sm %in% c("both", "lm"))
      g_line <- g_line +
      geom_smooth(method = "lm", linetype = "dashed", se = FALSE, linewidth = 0.8,
                  aes(weight = 1 / y_se^2, fill = after_scale(lighten(color, 0.5))))

    g_line <- g_line +
      geom_point(aes(size = n)) +
      scale_size_continuous(range = point_size_range, guide = "none",
                            breaks = if (is.null(point_size_breaks)) waiver() else point_size_breaks) +
      geom_errorbar(aes(ymin = pmax(y_mean - 1.96 * y_se, 0), ymax = y_mean + 1.96 * y_se),
                    width = 0, linewidth = 0.8, inherit.aes = TRUE) +
      scale_color_manual(values = palette_lc) +
      scale_x_continuous(breaks = scales::breaks_extended(n = 4)) +
      scale_y_continuous(breaks = scales::breaks_extended(n = 4)) +
      labs(x = if (!is.null(x_labels) && var %in% names(x_labels)) x_labels[[var]] else var,
           y = y_label_plot, color = NULL) +
      theme_classic(base_size = 16) +
      theme(legend.position   = c(0.3, 0.85),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key        = element_rect(fill = "transparent", color = NA))

    if (!is.na(y_lim)) g_line <- g_line + coord_cartesian(ylim = c(0, y_lim))
    plot(g_line)
    ggsave(file.path(output_directory_fig, paste0("line_", var, file_tag, "_bs_agg.pdf")),
           g_line, width = 4, height = 4, dpi = 300)
    ggsave(file.path(output_directory_fig, paste0("line_", var, file_tag, "_bs_agg.png")),
           g_line, width = 4, height = 4, dpi = 300)
  }
}

plot_lines_pooled(bins_all, vars = scatter_vars, n_bins = 15,
                  point_size_range = c(1, 5),
                  point_size_breaks = c(1000, 10000, 100000),
                  x_labels = list(
                    "r95p"        = expression(R95p ~ (days ~ yr^-1)),
                    "wet_sev_z"   = expression(Delta ~ "Precip. intensity (" * sigma * ")"),
                    "tree_cover"  = "Tree cover (%)",
                    "insect_mort" = expression(Infestation ~ frequency)
                  ),
                  smooth_method = list(
                    wet_sev_z   = "both",   # loess + linear
                    tree_cover  = "loess",  # loess only
                    r95p        = "loess",  # loess only
                    insect_mort = "lm"      # linear only
                  ),
                  y_max_by_var = list(insect_mort = NA),   # free y-axis (no clipping)
                  discrete_vars = c("insect_mort", "insect_defo"))

#### 2D: ppt x wet_sev_z coloured by mean mortality rate (matched samples, first bootstrap) ####
{
  clip <- function(df, var)
    dplyr::between(df[[var]], quantile(df[[var]], cut_lo, na.rm = TRUE),
                   quantile(df[[var]], cut_hi, na.rm = TRUE))

  matched_all <- as.data.frame(rbindlist(all_matched))
  plot_df_2d  <- matched_all %>% filter(clip(., "ppt"), clip(., "wet_sev_z"))

  mort_max <- 5

  shared_fill <- scale_fill_distiller(
    name = NULL, palette = "Reds", direction = 1,
    limits = c(0, mort_max), breaks = c(0, mort_max), labels = c("0", "≥ 5"),
    oob = scales::squish
  )

  make_panel <- function(lc) {
    ggplot(plot_df_2d %>% filter(lc_group == lc),
           aes(x = ppt, y = wet_sev_z, z = mean_mort_rate)) +
      stat_summary_2d(fun = mean, bins = 15) +
      shared_fill +
      labs(x = expression(MAP ~ (mm ~ yr^{-1})),
           y = expression(Delta ~ "Precip. intensity (" * sigma * ")")) +
      scale_x_continuous(breaks = scales::breaks_extended(n = 4)) +
      scale_y_continuous(breaks = scales::breaks_extended(n = 4)) +
      theme_classic(base_size = 16) +
      theme(legend.position = "none")
  }

  for (lc in c("Wetland", "Non-wetland")) {
    g_panel <- make_panel(lc)
    fname   <- paste0("2d_ppt_wetSev_mortality_", gsub("-", "", lc))
    ggsave(file.path(output_directory_fig, paste0(fname, ".pdf")),
           g_panel, width = 4, height = 4, dpi = 300, device = cairo_pdf)
    ggsave(file.path(output_directory_fig, paste0(fname, ".png")), g_panel, width = 4, height = 4, dpi = 300)
    plot(g_panel)
  }

  # Extract and save legend independently (horizontal, title on top)
  g_leg <- make_panel("Non-wetland") +
    guides(fill = guide_colorbar(
      title = expression(Mortality ~ rate ~ (ha^{-1} ~ yr^{-1})), title.position = "top", title.hjust = 1,
      barwidth = unit(2, "cm"), barheight = unit(0.3, "cm"), direction = "horizontal",
      ticks = FALSE, draw.ulim = FALSE, draw.llim = FALSE, frame.colour = "black", frame.linewidth = 0.3
    )) +
    theme(legend.position   = "bottom",
          legend.title      = element_text(size = 14, vjust = 0.8),
          legend.text       = element_text(size = 12, hjust = 1),
          legend.background = element_rect(fill = "transparent", color = NA),
          plot.background   = element_rect(fill = "transparent", color = NA),
          panel.background  = element_rect(fill = "transparent", color = NA))
  leg      <- ggpubr::get_legend(g_leg)
  leg_plot <- ggpubr::as_ggplot(leg)
  ggsave(file.path(output_directory_fig, "2d_ppt_wetSev_mortality_legend.pdf"),
         leg_plot, width = 2.5, height = 1.2, dpi = 300, bg = "transparent", device = cairo_pdf)
  ggsave(file.path(output_directory_fig, "2d_ppt_wetSev_mortality_legend.png"),
         leg_plot, width = 2.5, height = 1.2, dpi = 300, bg = "transparent")
  plot(leg_plot)
}
