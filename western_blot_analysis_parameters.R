################################################################################
# Western Blot Band & Lane Analysis in R
# ----------------------------------------
# This script:
#   1. Loads a Western blot image
#   2. Identifies horizontal bands (e.g., actin, tau) and saves a validation image
#   3. Divides each band into individual lanes and saves a validation image
#   4. Quantifies intensity for each band×lane ROI
#   5. Exports results to CSV
#
# REQUIRED PACKAGES: EBImage, magick, ggplot2, dplyr
#
# Install if needed:
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("EBImage")
#   install.packages(c("magick", "ggplot2", "dplyr"))
################################################################################

library(EBImage)
library(magick)
library(ggplot2)
library(dplyr)

# =============================================================================
# USER CONFIGURATION — EDIT THIS SECTION
# =============================================================================

# Path to your Western blot image (TIFF, PNG, JPEG, etc.)
IMAGE_PATH <- "my_western_blot.tif"

# Output directory (will be created if it doesn't exist)
OUTPUT_DIR <- "wb_analysis_output"

# --- Band definitions (horizontal rows in the blot) ---
# Define each band as a named list with y_start and y_end (pixel rows).
# TIP: Open the image in ImageJ or any viewer to find approximate pixel coordinates.
# The script will save a grayscale profile to help you refine these.
BANDS <- list(
  tau   = list(y_start = 80,  y_end = 160),
  actin = list(y_start = 250, y_end = 310)
)

# --- Lane definitions (vertical columns in the blot) ---
# Define lane boundaries as pixel x-coordinates (left edges), plus the final right edge.
# For example, 6 lanes would need 7 boundary values.
# You can also provide lane labels (treatment names).
LANE_BOUNDARIES <- c(30, 110, 190, 270, 350, 430, 510)

LANE_LABELS <- c(

  "Ctrl_1", "Ctrl_2", "Ctrl_3",
  "Treat_1", "Treat_2", "Treat_3"
)

# --- Analysis parameters ---
# Background subtraction method: "rolling_ball" or "margin"
# "margin" uses the edges of each ROI to estimate local background
# "rolling_ball" uses a simple local median filter approximation
BG_METHOD <- "margin"

# Margin width (pixels) for background estimation when using "margin" method
BG_MARGIN <- 5

# =============================================================================
# STEP 0: Setup — Create output directory, load and inspect image
# =============================================================================

cat("=" , rep("=", 70), "\n", sep = "")
cat("WESTERN BLOT ANALYSIS\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat("[OK] Created output directory:", OUTPUT_DIR, "\n")
}

# Load image
cat("[STEP 0] Loading image:", IMAGE_PATH, "\n")
if (!file.exists(IMAGE_PATH)) {
  stop("ERROR: Image file not found at '", IMAGE_PATH, "'. Please check the path.")
}

img <- readImage(IMAGE_PATH)
cat("  Image dimensions:", dim(img)[1], "x", dim(img)[2], "\n")
cat("  Color mode:", ifelse(colorMode(img) == 0, "Grayscale", "Color"), "\n")
cat("  Number of channels:", numberOfFrames(img, type = "render"), "\n")

# Convert to grayscale if needed
if (colorMode(img) != 0) {
  img_gray <- channel(img, "gray")
  cat("  Converted to grayscale.\n")
} else {
  img_gray <- img
}

# Save the raw grayscale image for reference
writeImage(img_gray, file.path(OUTPUT_DIR, "00_grayscale_input.png"))
cat("[SAVED] 00_grayscale_input.png\n\n")

# =============================================================================
# STEP 1: Generate intensity profiles to help define band/lane positions
# =============================================================================

cat("[STEP 1] Generating intensity profiles for band/lane positioning...\n")

img_matrix <- as.numeric(img_gray)
dim(img_matrix) <- dim(img_gray)[1:2]

# --- Horizontal profile (sum across columns → helps find BAND y-positions) ---
horizontal_profile <- rowMeans(img_matrix)

df_hprofile <- data.frame(
  y_pixel = seq_along(horizontal_profile),
  mean_intensity = horizontal_profile
)

p_hprofile <- ggplot(df_hprofile, aes(x = y_pixel, y = mean_intensity)) +
  geom_line(color = "steelblue", linewidth = 0.5) +
  # Overlay band regions
  annotate("rect",
    xmin = sapply(BANDS, `[[`, "y_start"),
    xmax = sapply(BANDS, `[[`, "y_end"),
    ymin = -Inf, ymax = Inf,
    fill = "red", alpha = 0.15
  ) +
  annotate("text",
    x = sapply(BANDS, function(b) mean(c(b$y_start, b$y_end))),
    y = max(horizontal_profile) * 0.95,
    label = names(BANDS),
    color = "red", fontface = "bold", size = 3.5
  ) +
  labs(
    title = "Horizontal Intensity Profile (rows)",
    subtitle = "Shaded regions = your defined bands. Peaks should align with bands.",
    x = "Y pixel (top → bottom)",
    y = "Mean intensity"
  ) +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "01_horizontal_profile.png"), p_hprofile,
  width = 10, height = 4, dpi = 150
)
cat("[SAVED] 01_horizontal_profile.png — Check that band regions align with peaks.\n")

# --- Vertical profile (sum across rows → helps find LANE x-positions) ---
vertical_profile <- colMeans(img_matrix)

df_vprofile <- data.frame(
  x_pixel = seq_along(vertical_profile),
  mean_intensity = vertical_profile
)

p_vprofile <- ggplot(df_vprofile, aes(x = x_pixel, y = mean_intensity)) +
  geom_line(color = "steelblue", linewidth = 0.5) +
  geom_vline(xintercept = LANE_BOUNDARIES, color = "red", linetype = "dashed", linewidth = 0.4) +
  labs(
    title = "Vertical Intensity Profile (columns)",
    subtitle = "Dashed red lines = your lane boundaries. Should fall between lane peaks.",
    x = "X pixel (left → right)",
    y = "Mean intensity"
  ) +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "01_vertical_profile.png"), p_vprofile,
  width = 10, height = 4, dpi = 150
)
cat("[SAVED] 01_vertical_profile.png — Check that lane boundaries are correct.\n\n")

# =============================================================================
# STEP 2: Draw rectangles around each BAND and save validation image
# =============================================================================

cat("[STEP 2] Drawing band-level ROIs on the image...\n")

# Use magick for annotation (much easier to draw colored rectangles with labels)
img_magick <- image_read(file.path(OUTPUT_DIR, "00_grayscale_input.png"))
img_info <- image_info(img_magick)
img_w <- img_info$width
img_h <- img_info$height

# Color palette for bands
band_colors <- c("red", "cyan", "yellow", "magenta", "green", "orange")

img_bands <- img_magick

for (i in seq_along(BANDS)) {
  band_name <- names(BANDS)[i]
  b <- BANDS[[i]]
  col <- band_colors[((i - 1) %% length(band_colors)) + 1]

  # Draw rectangle spanning full width around the band
  # magick uses geometry: widthxheight+x_offset+y_offset
  rect_x <- 0
  rect_y <- b$y_start
  rect_w <- img_w
  rect_h <- b$y_end - b$y_start

  # Create a transparent overlay with the rectangle
  draw_spec <- sprintf(
    "rectangle %d,%d %d,%d",
    rect_x, rect_y, rect_x + rect_w, rect_y + rect_h
  )

  img_bands <- image_draw(img_bands)
  rect(
    xleft = rect_x, ybottom = img_h - rect_y,
    xright = rect_x + rect_w, ytop = img_h - (rect_y + rect_h),
    border = col, lwd = 3, col = NA
  )
  # Add label
  text(
    x = 10, y = img_h - rect_y + 15,
    labels = toupper(band_name),
    col = col, cex = 1.5, font = 2, adj = c(0, 1)
  )
  dev.off()
}

image_write(img_bands, file.path(OUTPUT_DIR, "02_band_ROIs.png"))
cat("[SAVED] 02_band_ROIs.png — Rectangles around each band.\n\n")

# =============================================================================
# STEP 3: Draw rectangles for each LANE within each BAND
# =============================================================================

cat("[STEP 3] Drawing lane-level ROIs within each band...\n")

# Validate lane definitions
n_lanes <- length(LANE_BOUNDARIES) - 1
if (length(LANE_LABELS) != n_lanes) {
  stop(
    "ERROR: Number of LANE_LABELS (", length(LANE_LABELS),
    ") must equal number of lanes from LANE_BOUNDARIES (",
    n_lanes, ")."
  )
}

cat("  Number of lanes:", n_lanes, "\n")
cat("  Lane labels:", paste(LANE_LABELS, collapse = ", "), "\n")

img_lanes <- img_magick

# Store all ROI definitions for quantification
roi_list <- list()

img_lanes <- image_draw(img_lanes)

for (i in seq_along(BANDS)) {
  band_name <- names(BANDS)[i]
  b <- BANDS[[i]]
  col <- band_colors[((i - 1) %% length(band_colors)) + 1]

  for (j in seq_len(n_lanes)) {
    x_left  <- LANE_BOUNDARIES[j]
    x_right <- LANE_BOUNDARIES[j + 1]
    y_top   <- b$y_start
    y_bot   <- b$y_end

    # Draw rectangle (magick/R graphics coordinate: origin at top-left)
    rect(
      xleft = x_left, ybottom = img_h - y_top,
      xright = x_right, ytop = img_h - y_bot,
      border = col, lwd = 2, col = NA
    )

    # Small lane label inside the box
    text(
      x = (x_left + x_right) / 2,
      y = img_h - ((y_top + y_bot) / 2),
      labels = LANE_LABELS[j],
      col = col, cex = 0.7, font = 2
    )

    # Store ROI info
    roi_list[[length(roi_list) + 1]] <- data.frame(
      band       = band_name,
      lane       = LANE_LABELS[j],
      lane_index = j,
      x_left     = x_left,
      x_right    = x_right,
      y_top      = y_top,
      y_bottom   = y_bot,
      stringsAsFactors = FALSE
    )
  }

  # Band label on the left margin
  text(
    x = 5, y = img_h - b$y_start + 12,
    labels = toupper(band_name),
    col = col, cex = 1.2, font = 2, adj = c(0, 1)
  )
}

dev.off()

image_write(img_lanes, file.path(OUTPUT_DIR, "03_lane_ROIs.png"))
cat("[SAVED] 03_lane_ROIs.png — Individual lane rectangles within each band.\n\n")

roi_df <- bind_rows(roi_list)

# =============================================================================
# STEP 4: Quantify intensity in each ROI
# =============================================================================

cat("[STEP 4] Quantifying intensity for each band × lane ROI...\n")

# Helper: estimate background from margin pixels of an ROI
estimate_background_margin <- function(roi_pixels, margin) {
  nr <- nrow(roi_pixels)
  nc <- ncol(roi_pixels)
  if (nr <= 2 * margin || nc <= 2 * margin) {
    return(median(roi_pixels, na.rm = TRUE))
  }
  # Take pixels from the edges only
  top    <- roi_pixels[1:margin, ]
  bottom <- roi_pixels[(nr - margin + 1):nr, ]
  left   <- roi_pixels[, 1:margin]
  right  <- roi_pixels[, (nc - margin + 1):nc]
  margin_pixels <- c(as.numeric(top), as.numeric(bottom),
                     as.numeric(left), as.numeric(right))
  return(median(margin_pixels, na.rm = TRUE))
}

# Quantification results
results <- data.frame(
  band              = character(),
  lane              = character(),
  lane_index        = integer(),
  roi_width_px      = integer(),
  roi_height_px     = integer(),
  roi_area_px       = integer(),
  raw_integrated    = numeric(),
  raw_mean          = numeric(),
  bg_estimate       = numeric(),
  corrected_integrated = numeric(),
  corrected_mean    = numeric(),
  stringsAsFactors  = FALSE
)

for (r in seq_len(nrow(roi_df))) {
  roi <- roi_df[r, ]

  # Extract pixel region (note: EBImage matrix is [x, y])
  # Clamp to image bounds
  x1 <- max(1, roi$x_left + 1)
  x2 <- min(dim(img_gray)[1], roi$x_right)
  y1 <- max(1, roi$y_top + 1)
  y2 <- min(dim(img_gray)[2], roi$y_bottom)

  roi_pixels <- img_matrix[x1:x2, y1:y2]

  roi_w <- ncol(roi_pixels)
  roi_h <- nrow(roi_pixels)
  roi_area <- roi_w * roi_h

  raw_sum  <- sum(roi_pixels, na.rm = TRUE)
  raw_mean <- mean(roi_pixels, na.rm = TRUE)

  # Background estimation
  if (BG_METHOD == "margin") {
    bg <- estimate_background_margin(roi_pixels, BG_MARGIN)
  } else {
    # Simple rolling ball approximation: local median
    bg <- median(roi_pixels, na.rm = TRUE) * 0.5
  }

  corrected_pixels <- pmax(roi_pixels - bg, 0)
  corrected_sum  <- sum(corrected_pixels, na.rm = TRUE)
  corrected_mean <- mean(corrected_pixels, na.rm = TRUE)

  results <- rbind(results, data.frame(
    band              = roi$band,
    lane              = roi$lane,
    lane_index        = roi$lane_index,
    roi_width_px      = roi_w,
    roi_height_px     = roi_h,
    roi_area_px       = roi_area,
    raw_integrated    = round(raw_sum, 4),
    raw_mean          = round(raw_mean, 6),
    bg_estimate       = round(bg, 6),
    corrected_integrated = round(corrected_sum, 4),
    corrected_mean    = round(corrected_mean, 6),
    stringsAsFactors  = FALSE
  ))
}

cat("  Quantification complete for", nrow(results), "ROIs.\n\n")

# =============================================================================
# STEP 5: Save results to CSV
# =============================================================================

csv_path <- file.path(OUTPUT_DIR, "04_quantification_results.csv")
write.csv(results, csv_path, row.names = FALSE)
cat("[SAVED] 04_quantification_results.csv\n")

# Also print a summary table
cat("\n--- RESULTS SUMMARY ---\n\n")
print(results[, c("band", "lane", "raw_integrated", "corrected_integrated", "bg_estimate")],
  row.names = FALSE
)

# =============================================================================
# STEP 6: Generate a summary bar plot (optional but useful)
# =============================================================================

cat("\n[STEP 6] Generating summary bar plot...\n")

p_bars <- ggplot(results, aes(x = lane, y = corrected_integrated, fill = band)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  facet_wrap(~band, scales = "free_y") +
  labs(
    title = "Western Blot Quantification",
    subtitle = paste("Background method:", BG_METHOD),
    x = "Lane",
    y = "Background-corrected integrated intensity"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(OUTPUT_DIR, "05_quantification_barplot.png"), p_bars,
  width = 8, height = 5, dpi = 150
)
cat("[SAVED] 05_quantification_barplot.png\n")

# =============================================================================
# STEP 7: Normalized ratios (e.g., tau / actin per lane)
# =============================================================================

cat("\n[STEP 7] Computing normalized ratios...\n")

# Pivot to wide format: one row per lane, columns for each band's intensity
wide <- results %>%
  select(band, lane, lane_index, corrected_integrated) %>%
  tidyr::pivot_wider(
    names_from = band,
    values_from = corrected_integrated,
    names_prefix = "intensity_"
  )

# If both tau and actin exist, compute ratio
band_names_lower <- tolower(names(BANDS))
if ("tau" %in% band_names_lower && "actin" %in% band_names_lower) {
  wide <- wide %>%
    mutate(
      ratio_tau_over_actin = intensity_tau / intensity_actin
    )
  cat("  Computed tau/actin ratio per lane.\n")
}

ratio_csv <- file.path(OUTPUT_DIR, "06_normalized_ratios.csv")
write.csv(wide, ratio_csv, row.names = FALSE)
cat("[SAVED] 06_normalized_ratios.csv\n")

cat("\n--- NORMALIZED DATA ---\n\n")
print(as.data.frame(wide), row.names = FALSE)

# =============================================================================
# DONE
# =============================================================================

cat("\n", rep("=", 72), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 72), "\n", sep = "")
cat("\nAll outputs saved to:", normalizePath(OUTPUT_DIR), "\n\n")
cat("Output files:\n")
cat("  00_grayscale_input.png        — Grayscale version of your input\n")
cat("  01_horizontal_profile.png     — Y-axis intensity profile (validate band positions)\n")
cat("  01_vertical_profile.png       — X-axis intensity profile (validate lane positions)\n")
cat("  02_band_ROIs.png              — Image with rectangles around each band\n")
cat("  03_lane_ROIs.png              — Image with rectangles for each lane × band\n")
cat("  04_quantification_results.csv — Raw & corrected intensities per ROI\n")
cat("  05_quantification_barplot.png — Bar plot of intensities\n")
cat("  06_normalized_ratios.csv      — Ratios (e.g., tau/actin) per lane\n")
cat("\nTIP: If band/lane positions look wrong in the validation images,\n")
cat("     adjust BANDS and LANE_BOUNDARIES in the config section and re-run.\n")
