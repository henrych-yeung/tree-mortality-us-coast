# ALIGN TEMPORAL FIGURES — stack the group-comparison figure (3a) above the
# driver overview figure (3c), with matched facet/column widths.
# Run AFTER 3a_lowElev_temporal.R and 3c_plot_temporal_driver.R (both save
# their plot objects as .rds into the shared temporal output folder).
# Author: Henry CH Yeung, UVA

library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)

base_dir <- "/Volumes/Henry2/ghostForest_us/_publish_product"   # EDIT to your local path
setwd(base_dir)

CONF_THRESH_STR <- "conf0d5"   # must match the value used in 3a & 3c

fig_dir <- paste0("_output/3_temporal_", CONF_THRESH_STR)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

group_tag <- "none"   # match bar_mode in 3a_lowElev_temporal.R ("area" | "hotspot" | "none")
g_group  <- readRDS(file.path(fig_dir, paste0("g_group_", group_tag, ".rds")))
g_driver <- readRDS(file.path(fig_dir, "g_driver.rds"))

# Drop x-axis title from top panel (shared with bottom)
g_group <- g_group + theme(axis.title.x = element_blank())
# Drop facet strip labels from bottom panel (already shown on top)
g_driver <- g_driver + theme(strip.text = element_blank())

# Standardize y-axis label width in g_group to match g_driver.
# g_driver labels are like "0.20" (4 chars); g_group has "0"-"10" (1-2 chars).
# Right-padding integers to width 4 ("0   "-"10  ") keeps rendered label width
# consistent so unit.pmax on the gtable columns aligns correctly.
# Preserve the existing sec.axis (hotspot/area %) if present -- adding a new
# scale_y_continuous would otherwise wipe it out. For bar_mode == "none"
# there's no sec.axis to preserve.
y_scale <- Filter(function(s) "y" %in% s$aesthetics, g_group$scales$scales)
existing_sec <- if (length(y_scale)) y_scale[[1]]$secondary.axis else waiver()
g_group <- g_group +
  scale_y_continuous(
    labels   = function(x) formatC(as.integer(x), width = 4, flag = "-"),
    breaks   = scales::pretty_breaks(n = 4),
    sec.axis = existing_sec
  )

# Convert to gtable and equalize column widths so panels line up
gt_group  <- ggplotGrob(g_group)
gt_driver <- ggplotGrob(g_driver)
shared_w  <- unit.pmax(gt_group$widths, gt_driver$widths)
gt_group$widths  <- shared_w
gt_driver$widths <- shared_w

white_spacer <- rectGrob(gp = gpar(col = NA, fill = "white"))

g_stacked <- arrangeGrob(
  gt_group, white_spacer, gt_driver,
  ncol = 1,
  heights = unit(c(1, 0.1, 1), "null")
)

ggsave(file.path(fig_dir, "temporal_aligned.pdf"), g_stacked, width = 16, height = 8.2, dpi = 300)
ggsave(file.path(fig_dir, "temporal_aligned.png"), g_stacked, width = 16, height = 8.2, dpi = 300)
grid.newpage(); grid.draw(g_stacked)
cat("Saved to:", fig_dir, "\n")
