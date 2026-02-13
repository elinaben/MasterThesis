################################################################################
# Automatic Western Blot Band & Lane Detection
# -----------------------------------------------
# Fully automated — no pixel coordinates needed.
# Detects bands and lanes from intensity profiles, marks them,
# and outputs a matrix of intensities to CSV.
#
# REQUIRED PACKAGES: EBImage, magick, ggplot2, dplyr, tidyr
#
# Install if needed:
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("EBImage")
#   install.packages(c("magick", "ggplot2", "dplyr", "tidyr"))
################################################################################

library(EBImage)
library(magick)
library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# USER CONFIGURATION — minimal settings
# =============================================================================

# Path to your Western blot image
IMAGE_PATH <- "figs/ELINA_2026-02-06_13h04m41s_IRDye_680RD_.jpg"

# Output directory
OUTPUT_DIR <- "wb_auto_output"

# --- Tuning knobs (sensible defaults — adjust only if detection is off) ---

# Should the image be inverted? (TRUE if bands are DARK on LIGHT background,
# which is the most common case for chemiluminescent / Coomassie blots)
INVERT_IMAGE <- TRUE

# Smoothing window for intensity profiles (pixels). Larger = less noise, but
# may merge close bands/lanes. Smaller = more sensitive but noisier.
SMOOTH_WINDOW_BANDS <- 15
SMOOTH_WINDOW_LANES <- 15

# Peak prominence threshold as a fraction of the profile range (0–1).
# Higher = only strong bands/lanes detected. Lower = more sensitive.
PEAK_THRESHOLD_BANDS <- 0.15
PEAK_THRESHOLD_LANES <- 0.10

# Minimum width (pixels) for a detected band or lane to be kept
MIN_BAND_HEIGHT <- 10
MIN_LANE_WIDTH  <- 10

# Padding (pixels) to add around detected peak centers to define ROI edges
# Set to 0 to use the valley-to-valley boundaries instead (recommended)
PADDING <- 0

# Crop margins — fraction of image to ignore at edges (0–0.1 range)
# Useful to exclude labels or artifacts at the borders
CROP_TOP    <- 0.0
CROP_BOTTOM <- 0.0
CROP_LEFT   <- 0.0
CROP_RIGHT  <- 0.0

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Simple moving average smoother
smooth_profile <- function(x, window) {
  if (window <= 1) return(x)
  k <- floor(window / 2)
  n <- length(x)
  out <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1, i - k)
    hi <- min(n, i + k)
    out[i] <- mean(x[lo:hi])
  }
  out
}

#' Detect peaks in a 1D profile.
#' Returns a data.frame with columns: center, height, left, right
#' where left/right are the valley boundaries on each side.
detect_peaks <- function(profile, threshold_frac = 0.15, min_width = 10,
                         smooth_window = 15) {
  
  prof <- smooth_profile(profile, smooth_window)
  n <- length(prof)
  
  # Threshold: only consider peaks above this
  prof_range <- max(prof) - min(prof)
  threshold  <- min(prof) + prof_range * threshold_frac
  
  # Find local maxima
  is_peak <- rep(FALSE, n)
  for (i in 2:(n - 1)) {
    if (prof[i] > prof[i - 1] && prof[i] > prof[i + 1] && prof[i] >= threshold) {
      is_peak[i] <- TRUE
    }
  }
  
  peak_positions <- which(is_peak)
  if (length(peak_positions) == 0) return(NULL)
  
  # For each peak, find the valley boundaries (descend left and right)
  results <- list()
  
  for (pk in peak_positions) {
    # Walk left to find valley
    left_bound <- pk
    for (j in (pk - 1):1) {
      if (prof[j] > prof[j + 1]) {
        left_bound <- j + 1
        break
      }
      left_bound <- j
    }
    
    # Walk right to find valley
    right_bound <- pk
    for (j in (pk + 1):n) {
      if (prof[j] > prof[j - 1]) {
        right_bound <- j - 1
        break
      }
      right_bound <- j
    }
    
    width <- right_bound - left_bound + 1
    if (width >= min_width) {
      results[[length(results) + 1]] <- data.frame(
        center = pk,
        height = prof[pk],
        left   = left_bound,
        right  = right_bound,
        width  = width
      )
    }
  }
  
  if (length(results) == 0) return(NULL)
  peaks_df <- bind_rows(results)
  
  # Merge overlapping peaks (keep the taller one's center, expand boundaries)
  if (nrow(peaks_df) > 1) {
    peaks_df <- peaks_df[order(peaks_df$left), ]
    merged <- list(peaks_df[1, ])
    for (i in 2:nrow(peaks_df)) {
      prev <- merged[[length(merged)]]
      curr <- peaks_df[i, ]
      if (curr$left <= prev$right) {
        # Overlapping — merge
        new_left  <- min(prev$left, curr$left)
        new_right <- max(prev$right, curr$right)
        winner    <- if (curr$height > prev$height) curr else prev
        merged[[length(merged)]] <- data.frame(
          center = winner$center,
          height = max(prev$height, curr$height),
          left   = new_left,
          right  = new_right,
          width  = new_right - new_left + 1
        )
      } else {
        merged[[length(merged) + 1]] <- curr
      }
    }
    peaks_df <- bind_rows(merged)
  }
  
  peaks_df
}

# =============================================================================
# STEP 0: Setup and load image
# =============================================================================

cat(strrep("=", 72), "\n")
cat("AUTOMATIC WESTERN BLOT ANALYSIS\n")
cat(strrep("=", 72), "\n\n")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("[STEP 0] Loading image:", IMAGE_PATH, "\n")
if (!file.exists(IMAGE_PATH)) {
  stop("Image not found: '", IMAGE_PATH, "'")
}

img_raw <- readImage(IMAGE_PATH)
cat("  Raw dimensions: ", dim(img_raw)[1], " x ", dim(img_raw)[2], "\n", sep = "")

# Convert to grayscale
if (colorMode(img_raw) != 0) {
  img_gray <- channel(img_raw, "gray")
} else {
  img_gray <- img_raw
}

# Invert if needed (so bands become bright on dark background)
if (INVERT_IMAGE) {
  img_gray <- 1 - img_gray
  cat("  Image inverted (bands are now bright).\n")
}

# Crop margins
nr <- dim(img_gray)[1]  # Note: EBImage stores as [x, y] i.e. [width, height]
nc <- dim(img_gray)[2]
x1 <- max(1, round(nr * CROP_LEFT) + 1)
x2 <- min(nr, round(nr * (1 - CROP_RIGHT)))
y1 <- max(1, round(nc * CROP_TOP) + 1)
y2 <- min(nc, round(nc * (1 - CROP_BOTTOM)))
img_gray <- img_gray[x1:x2, y1:y2]

cat("  Working dimensions: ", dim(img_gray)[1], " x ", dim(img_gray)[2], "\n", sep = "")

writeImage(normalize(img_gray), file.path(OUTPUT_DIR, "00_preprocessed.png"))
cat("[SAVED] 00_preprocessed.png\n\n")

# Get raw matrix (columns = x, rows = y for standard image orientation)
# EBImage: dim1 = x (width), dim2 = y (height)
img_mat <- as.numeric(img_gray)
dim(img_mat) <- dim(img_gray)[1:2]
img_w <- dim(img_mat)[1]
img_h <- dim(img_mat)[2]

# =============================================================================
# STEP 1: Detect BANDS (horizontal rows)
# =============================================================================

cat("[STEP 1] Detecting bands (horizontal signal rows)...\n")

# Average intensity per row (along y-axis = dim2)
# For each y position, average across all x positions
band_profile <- colMeans(img_mat)  # length = img_h

bands <- detect_peaks(
  profile       = band_profile,
  threshold_frac = PEAK_THRESHOLD_BANDS,
  min_width     = MIN_BAND_HEIGHT,
  smooth_window = SMOOTH_WINDOW_BANDS
)

if (is.null(bands) || nrow(bands) == 0) {
  stop("No bands detected! Try lowering PEAK_THRESHOLD_BANDS or SMOOTH_WINDOW_BANDS.")
}

bands <- bands[order(bands$center), ]
bands$band_id <- seq_len(nrow(bands))
bands$band_name <- paste0("Band_", bands$band_id)

cat("  Detected", nrow(bands), "bands:\n")
for (i in seq_len(nrow(bands))) {
  cat(sprintf("    %s: y = %d–%d (center: %d, width: %d px)\n",
              bands$band_name[i], bands$left[i], bands$right[i],
              bands$center[i], bands$width[i]
  ))
}

# Save band profile plot
df_bp <- data.frame(y = seq_along(band_profile), intensity = band_profile)
smooth_bp <- smooth_profile(band_profile, SMOOTH_WINDOW_BANDS)
df_bp$smoothed <- smooth_bp

p_bands <- ggplot(df_bp, aes(x = y)) +
  geom_line(aes(y = intensity), color = "gray60", linewidth = 0.3) +
  geom_line(aes(y = smoothed), color = "steelblue", linewidth = 0.8) +
  geom_vline(xintercept = bands$center, color = "red", linetype = "solid", linewidth = 0.5) +
  annotate("rect",
           xmin = bands$left, xmax = bands$right,
           ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.12
  ) +
  annotate("text",
           x = bands$center, y = max(df_bp$smoothed) * 1.02,
           label = bands$band_name, color = "red", size = 3, angle = 0, vjust = 0
  ) +
  labs(
    title = "Band Detection — Horizontal Intensity Profile",
    subtitle = sprintf("Detected %d band(s). Blue = smoothed profile. Red = detected bands.", nrow(bands)),
    x = "Y position (pixels, top → bottom)",
    y = "Mean intensity"
  ) +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "01_band_detection_profile.png"), p_bands,
       width = 10, height = 4, dpi = 150
)
cat("[SAVED] 01_band_detection_profile.png\n\n")

# =============================================================================
# STEP 2: Detect LANES (vertical columns)
# =============================================================================

cat("[STEP 2] Detecting lanes (vertical signal columns)...\n")

# Average intensity per column (along x-axis = dim1)
# Only use pixels within detected band regions for better lane detection
band_mask <- rep(FALSE, img_h)
for (i in seq_len(nrow(bands))) {
  band_mask[bands$left[i]:bands$right[i]] <- TRUE
}

if (sum(band_mask) > 0) {
  lane_profile <- rowMeans(img_mat[, band_mask, drop = FALSE])
} else {
  lane_profile <- rowMeans(img_mat)
}

lanes <- detect_peaks(
  profile       = lane_profile,
  threshold_frac = PEAK_THRESHOLD_LANES,
  min_width     = MIN_LANE_WIDTH,
  smooth_window = SMOOTH_WINDOW_LANES
)

if (is.null(lanes) || nrow(lanes) == 0) {
  stop("No lanes detected! Try lowering PEAK_THRESHOLD_LANES or SMOOTH_WINDOW_LANES.")
}

lanes <- lanes[order(lanes$center), ]
lanes$lane_id <- seq_len(nrow(lanes))
lanes$lane_name <- paste0("Lane_", lanes$lane_id)

cat("  Detected", nrow(lanes), "lanes:\n")
for (i in seq_len(nrow(lanes))) {
  cat(sprintf("    %s: x = %d–%d (center: %d, width: %d px)\n",
              lanes$lane_name[i], lanes$left[i], lanes$right[i],
              lanes$center[i], lanes$width[i]
  ))
}

# Save lane profile plot
df_lp <- data.frame(x = seq_along(lane_profile), intensity = lane_profile)
smooth_lp <- smooth_profile(lane_profile, SMOOTH_WINDOW_LANES)
df_lp$smoothed <- smooth_lp

p_lanes <- ggplot(df_lp, aes(x = x)) +
  geom_line(aes(y = intensity), color = "gray60", linewidth = 0.3) +
  geom_line(aes(y = smoothed), color = "steelblue", linewidth = 0.8) +
  geom_vline(xintercept = lanes$center, color = "darkgreen", linetype = "solid", linewidth = 0.5) +
  annotate("rect",
           xmin = lanes$left, xmax = lanes$right,
           ymin = -Inf, ymax = Inf,
           fill = "green", alpha = 0.12
  ) +
  annotate("text",
           x = lanes$center, y = max(df_lp$smoothed) * 1.02,
           label = lanes$lane_name, color = "darkgreen", size = 3, angle = 45, vjust = 0, hjust = 0
  ) +
  labs(
    title = "Lane Detection — Vertical Intensity Profile",
    subtitle = sprintf("Detected %d lane(s). Green = detected lanes. Profile computed within band regions.", nrow(lanes)),
    x = "X position (pixels, left → right)",
    y = "Mean intensity"
  ) +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "02_lane_detection_profile.png"), p_lanes,
       width = 10, height = 4, dpi = 150
)
cat("[SAVED] 02_lane_detection_profile.png\n\n")

# =============================================================================
# STEP 3: Validation image — rectangles around each detected BAND
# =============================================================================

cat("[STEP 3] Drawing detected band ROIs on image...\n")

img_for_draw <- normalize(img_gray)
writeImage(img_for_draw, file.path(OUTPUT_DIR, "tmp_draw.png"))
img_magick <- image_read(file.path(OUTPUT_DIR, "tmp_draw.png"))
info <- image_info(img_magick)
draw_w <- info$width
draw_h <- info$height

band_colors <- c("#FF3333", "#33CCFF", "#FFCC00", "#FF66FF", "#33FF66", "#FF9933")

img_draw <- image_draw(img_magick)

for (i in seq_len(nrow(bands))) {
  col <- band_colors[((i - 1) %% length(band_colors)) + 1]
  y_top <- bands$left[i]
  y_bot <- bands$right[i]
  
  rect(
    xleft = 0, xright = draw_w,
    # image_draw uses standard graphics coords: origin top-left, y increases down
    ytop = y_top, ybottom = y_bot,
    border = col, lwd = 3, col = NA
  )
  text(
    x = 5, y = y_top - 3,
    labels = bands$band_name[i],
    col = col, cex = 1.3, font = 2, adj = c(0, 1)
  )
}
dev.off()

image_write(img_draw, file.path(OUTPUT_DIR, "03_detected_bands.png"))
cat("[SAVED] 03_detected_bands.png\n\n")

# =============================================================================
# STEP 4: Validation image — rectangles for each LANE × BAND cell
# =============================================================================

cat("[STEP 4] Drawing lane × band grid on image...\n")

img_draw2 <- image_draw(img_magick)

for (i in seq_len(nrow(bands))) {
  col <- band_colors[((i - 1) %% length(band_colors)) + 1]
  y_top <- bands$left[i]
  y_bot <- bands$right[i]
  
  for (j in seq_len(nrow(lanes))) {
    x_left  <- lanes$left[j]
    x_right <- lanes$right[j]
    
    rect(
      xleft = x_left, xright = x_right,
      ytop = y_top, ybottom = y_bot,
      border = col, lwd = 2, col = NA
    )
    
    # Label: B1L1, B1L2, etc.
    cell_label <- paste0("B", i, "L", j)
    text(
      x = (x_left + x_right) / 2,
      y = (y_top + y_bot) / 2,
      labels = cell_label,
      col = col, cex = 0.65, font = 2
    )
  }
  
  # Band label on the right side
  text(
    x = draw_w - 5, y = (y_top + y_bot) / 2,
    labels = bands$band_name[i],
    col = col, cex = 1.0, font = 2, adj = c(1, 0.5)
  )
}

# Lane labels at the top
for (j in seq_len(nrow(lanes))) {
  text(
    x = (lanes$left[j] + lanes$right[j]) / 2,
    y = max(0, min(bands$left) - 15),
    labels = lanes$lane_name[j],
    col = "white", cex = 0.75, font = 2
  )
}

dev.off()

image_write(img_draw2, file.path(OUTPUT_DIR, "04_detected_grid.png"))
cat("[SAVED] 04_detected_grid.png\n\n")

# =============================================================================
# STEP 5: Quantify intensity for every cell in the grid
# =============================================================================

cat("[STEP 5] Quantifying intensity for each band × lane cell...\n")

results <- data.frame(
  band_id    = integer(),
  lane_id    = integer(),
  cell_label = character(),
  x_left     = integer(),
  x_right    = integer(),
  y_top      = integer(),
  y_bottom   = integer(),
  area_px    = integer(),
  raw_integrated_intensity    = numeric(),
  raw_mean_intensity          = numeric(),
  bg_median                   = numeric(),
  corrected_integrated_intensity = numeric(),
  corrected_mean_intensity    = numeric(),
  max_pixel_intensity         = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(bands))) {
  y1 <- bands$left[i]
  y2 <- bands$right[i]
  
  for (j in seq_len(nrow(lanes))) {
    x1 <- lanes$left[j]
    x2 <- lanes$right[j]
    
    # Extract ROI pixels: img_mat[x, y]
    roi <- img_mat[
      max(1, x1):min(img_w, x2),
      max(1, y1):min(img_h, y2),
      drop = FALSE
    ]
    
    roi_area   <- length(roi)
    raw_sum    <- sum(roi)
    raw_mean   <- mean(roi)
    max_px     <- max(roi)
    
    # Background: median of the border pixels (2-pixel margin)
    nr <- nrow(roi)
    nc <- ncol(roi)
    margin_px <- c(
      as.numeric(roi[1:min(2, nr), ]),
      as.numeric(roi[max(1, nr - 1):nr, ]),
      as.numeric(roi[, 1:min(2, nc)]),
      as.numeric(roi[, max(1, nc - 1):nc])
    )
    bg <- median(margin_px)
    
    corrected_roi  <- pmax(roi - bg, 0)
    corrected_sum  <- sum(corrected_roi)
    corrected_mean <- mean(corrected_roi)
    
    results <- rbind(results, data.frame(
      band_id    = i,
      lane_id    = j,
      cell_label = paste0("B", i, "L", j),
      x_left     = x1,
      x_right    = x2,
      y_top      = y1,
      y_bottom   = y2,
      area_px    = roi_area,
      raw_integrated_intensity    = round(raw_sum, 4),
      raw_mean_intensity          = round(raw_mean, 6),
      bg_median                   = round(bg, 6),
      corrected_integrated_intensity = round(corrected_sum, 4),
      corrected_mean_intensity    = round(corrected_mean, 6),
      max_pixel_intensity         = round(max_px, 6),
      stringsAsFactors = FALSE
    ))
  }
}

# Save detailed results
write.csv(results, file.path(OUTPUT_DIR, "05_detailed_results.csv"), row.names = FALSE)
cat("[SAVED] 05_detailed_results.csv\n")

# =============================================================================
# STEP 6: Output MATRIX format (bands as rows, lanes as columns)
# =============================================================================

cat("[STEP 6] Creating intensity matrix (bands × lanes)...\n")

# --- Raw integrated intensity matrix ---
mat_raw <- results %>%
  select(band_id, lane_id, raw_integrated_intensity) %>%
  pivot_wider(
    names_from  = lane_id,
    values_from = raw_integrated_intensity,
    names_prefix = "Lane_"
  ) %>%
  mutate(band_id = paste0("Band_", band_id)) %>%
  rename(Band = band_id)

write.csv(mat_raw, file.path(OUTPUT_DIR, "06_matrix_raw_integrated.csv"), row.names = FALSE)
cat("[SAVED] 06_matrix_raw_integrated.csv\n")

# --- Corrected integrated intensity matrix ---
mat_corr <- results %>%
  select(band_id, lane_id, corrected_integrated_intensity) %>%
  pivot_wider(
    names_from  = lane_id,
    values_from = corrected_integrated_intensity,
    names_prefix = "Lane_"
  ) %>%
  mutate(band_id = paste0("Band_", band_id)) %>%
  rename(Band = band_id)

write.csv(mat_corr, file.path(OUTPUT_DIR, "06_matrix_corrected_integrated.csv"), row.names = FALSE)
cat("[SAVED] 06_matrix_corrected_integrated.csv\n")

# --- Corrected mean intensity matrix ---
mat_mean <- results %>%
  select(band_id, lane_id, corrected_mean_intensity) %>%
  pivot_wider(
    names_from  = lane_id,
    values_from = corrected_mean_intensity,
    names_prefix = "Lane_"
  ) %>%
  mutate(band_id = paste0("Band_", band_id)) %>%
  rename(Band = band_id)

write.csv(mat_mean, file.path(OUTPUT_DIR, "06_matrix_corrected_mean.csv"), row.names = FALSE)
cat("[SAVED] 06_matrix_corrected_mean.csv\n\n")

# Print matrices to console
cat("--- RAW INTEGRATED INTENSITY MATRIX ---\n")
print(as.data.frame(mat_raw), row.names = FALSE)

cat("\n--- CORRECTED INTEGRATED INTENSITY MATRIX ---\n")
print(as.data.frame(mat_corr), row.names = FALSE)

cat("\n--- CORRECTED MEAN INTENSITY MATRIX ---\n")
print(as.data.frame(mat_mean), row.names = FALSE)

# =============================================================================
# STEP 7: Summary bar plot
# =============================================================================

cat("\n[STEP 7] Generating summary plots...\n")

results$band_label <- paste0("Band_", results$band_id)
results$lane_label <- paste0("Lane_", results$lane_id)

p_summary <- ggplot(results,
                    aes(x = lane_label, y = corrected_integrated_intensity, fill = band_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.3) +
  facet_wrap(~band_label, scales = "free_y", ncol = 1) +
  labs(
    title = "Automatic Western Blot Quantification",
    x = "Lane", y = "Background-corrected integrated intensity"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(OUTPUT_DIR, "07_summary_barplot.png"), p_summary,
       width = 8, height = 3 + 2.5 * nrow(bands), dpi = 150
)
cat("[SAVED] 07_summary_barplot.png\n")

# Heatmap of the matrix
mat_long <- results %>%
  select(band_label, lane_label, corrected_integrated_intensity)

p_heat <- ggplot(mat_long,
                 aes(x = lane_label, y = band_label, fill = corrected_integrated_intensity)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = round(corrected_integrated_intensity, 1)),
            color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Intensity Heatmap (Band × Lane)",
    x = "Lane", y = "Band", fill = "Corrected\nIntensity"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

ggsave(file.path(OUTPUT_DIR, "07_heatmap.png"), p_heat,
       width = max(6, 1.2 * nrow(lanes)), height = max(3, 1.5 * nrow(bands)), dpi = 150
)
cat("[SAVED] 07_heatmap.png\n")

# Cleanup temp file
file.remove(file.path(OUTPUT_DIR, "tmp_draw.png"))

# =============================================================================
# DONE
# =============================================================================

cat("\n", strrep("=", 72), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 72), "\n\n")
cat("Detected: ", nrow(bands), " bands × ", nrow(lanes), " lanes = ",
    nrow(bands) * nrow(lanes), " ROIs\n\n", sep = "")
cat("All outputs in:", normalizePath(OUTPUT_DIR), "\n\n")
cat("Files:\n")
cat("  00_preprocessed.png                 — Grayscale (optionally inverted) input\n")
cat("  01_band_detection_profile.png       — Horizontal profile + detected bands\n")
cat("  02_lane_detection_profile.png       — Vertical profile + detected lanes\n")
cat("  03_detected_bands.png               — Validation: band rectangles on image\n")
cat("  04_detected_grid.png                — Validation: full band × lane grid\n")
cat("  05_detailed_results.csv             — All metrics per ROI (long format)\n")
cat("  06_matrix_raw_integrated.csv        — Band × Lane matrix (raw intensity)\n")
cat("  06_matrix_corrected_integrated.csv  — Band × Lane matrix (BG-corrected)\n")
cat("  06_matrix_corrected_mean.csv        — Band × Lane matrix (corrected mean)\n")
cat("  07_summary_barplot.png              — Bar plot per band\n")
cat("  07_heatmap.png                      — Heatmap of the intensity matrix\n")
cat("\nTuning tips:\n")
cat("  - Too many bands/lanes? → Increase PEAK_THRESHOLD or SMOOTH_WINDOW\n")
cat("  - Missing bands/lanes?  → Decrease PEAK_THRESHOLD or SMOOTH_WINDOW\n")
cat("  - Merging neighbors?    → Decrease SMOOTH_WINDOW\n")
cat("  - Edge artifacts?       → Increase CROP_* margins\n")