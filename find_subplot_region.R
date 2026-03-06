# ─────────────────────────────────────────────────────
# Western Blot Lane Area Calculator — R version
# ─────────────────────────────────────────────────────
# Requirements:
#   install.packages(c("tiff", "pracma", "dplyr"))
#
# Usage (from R):
#   source("analyze_wb_plots.R")
#   results <- analyze_files(c("plot1.tif", "plot2.tif"), max_bands = 2)
#   write.csv(results, "western_blot_areas.csv", row.names = FALSE)
#
# Or run from the terminal:
#   Rscript analyze_wb_plots.R plot1.tif plot2.tif --bands 2
# ─────────────────────────────────────────────────────

library(tiff)     # readTIFF()
library(pracma)   # findpeaks()
library(dplyr)

# ── 1. Detect subplot (lane) row boundaries ───────────────────────────────────
find_subplot_regions <- function(arr, min_height = 20) {
  # arr is a matrix [row, col], values 0–1 (tiff package default)
  row_means  <- rowMeans(arr)
  white_rows <- row_means > (250 / 255)
  
  regions  <- list()
  in_plot  <- FALSE
  start    <- NULL
  
  for (i in seq_along(white_rows)) {
    if (!white_rows[i] && !in_plot) {
      in_plot <- TRUE
      start   <- i
    } else if (white_rows[i] && in_plot) {
      in_plot <- FALSE
      if ((i - 1 - start) > min_height)
        regions <- c(regions, list(c(start, i - 1)))
    }
  }
  # Handle subplot that runs to the bottom of the image
  if (in_plot && (nrow(arr) - start) > min_height)
    regions <- c(regions, list(c(start, nrow(arr))))
  
  regions
}

# ── 2. Trace the curve (find dark-pixel y-position per column) ────────────────
trace_curve <- function(subplot,
                        top_skip  = 5,  bottom_skip = 3,
                        side_skip = 3,  dark_threshold = 10 / 255) {
  h        <- nrow(subplot)
  interior <- subplot[(top_skip + 1):(h - bottom_skip),
                      (side_skip + 1):(ncol(subplot) - side_skip)]
  
  n_cols  <- ncol(interior)
  curve_y <- rep(NA_real_, n_cols)
  
  for (col in seq_len(n_cols)) {
    col_data  <- interior[, col]
    dark_rows <- which(col_data < dark_threshold)
    if (length(dark_rows) > 0)
      curve_y[col] <- mean(dark_rows)   # centre of the line
  }
  
  valid <- which(!is.na(curve_y))
  if (length(valid) == 0) return(NULL)
  
  list(col_idx = valid, y_vals = curve_y[valid])
}

# ── 3. Area between curve segment and its linear baseline ─────────────────────
area_under_segment <- function(y_vals) {
  n        <- length(y_vals)
  baseline <- seq(y_vals[1], y_vals[n], length.out = n)
  signal   <- pmax(y_vals - baseline, 0)   # only downward deflection counts
  # Trapezoidal integration (spacing = 1 pixel)
  sum(diff(seq_len(n)) * (head(signal, -1) + tail(signal, -1)) / 2)
}

# ── 4. Detect band boundaries (split at valley between peaks) ─────────────────
detect_bands <- function(y_vals, max_bands = 2) {
  n       <- length(y_vals)
  min_dist <- max(floor(n / 6), 5)
  
  # findpeaks() from pracma — returns matrix [peak_val, peak_idx, start, end]
  peaks_mat <- findpeaks(y_vals, minpeakdistance = min_dist, minpeakheight = 5)
  
  if (is.null(peaks_mat) || nrow(peaks_mat) < 2) {
    return(list(c(1L, n)))    # single band
  }
  
  # Keep the `max_bands` most prominent peaks (highest values)
  prominence_order <- order(peaks_mat[, 1], decreasing = TRUE)
  selected         <- peaks_mat[prominence_order[seq_len(min(max_bands, nrow(peaks_mat)))], , drop = FALSE]
  peak_cols        <- sort(selected[, 2])  # column positions, left-to-right
  
  # Split at valley minimum between each consecutive peak pair
  splits <- c(1L)
  for (i in seq_len(length(peak_cols) - 1)) {
    seg_start  <- peak_cols[i]
    seg_end    <- peak_cols[i + 1]
    valley_pos <- seg_start + which.min(y_vals[seg_start:seg_end]) - 1L
    splits     <- c(splits, valley_pos)
  }
  splits <- c(splits, n)
  
  # Return list of (start, end) pairs
  lapply(seq_len(length(splits) - 1), function(i) c(splits[i], splits[i + 1]))
}

# ── 5. Analyse one TIFF file ──────────────────────────────────────────────────
analyze_image <- function(tif_path, max_bands = 2) {
  arr      <- readTIFF(tif_path, as.is = FALSE)   # values in [0, 1]
  # If image is RGB, convert to grayscale (luminance)
  if (length(dim(arr)) == 3) arr <- arr[,,1] * 0.299 +
      arr[,,2] * 0.587 +
      arr[,,3] * 0.114
  
  regions <- find_subplot_regions(arr)
  message(sprintf("  %s: %d lanes detected", basename(tif_path), length(regions)))
  
  results <- vector("list", length(regions) * max_bands)
  row_n   <- 0L
  
  for (lane_idx in seq_along(regions)) {
    r       <- regions[[lane_idx]]
    subplot <- arr[r[1]:r[2], ]
    curve   <- trace_curve(subplot)
    
    if (is.null(curve)) {
      message(sprintf("  [WARN] Lane %d: no curve found, skipping", lane_idx))
      next
    }
    
    bands <- detect_bands(curve$y_vals, max_bands = max_bands)
    
    for (band_idx in seq_along(bands)) {
      seg     <- bands[[band_idx]]
      segment <- curve$y_vals[seg[1]:seg[2]]
      area    <- area_under_segment(segment)
      
      band_label <- if (length(bands) > 1) paste0("band_", band_idx) else "single"
      message(sprintf("  Lane %d | %s: area = %.2f px\u00b2", lane_idx, band_label, area))
      
      row_n <- row_n + 1L
      results[[row_n]] <- data.frame(
        image            = basename(tif_path),
        lane             = lane_idx,
        band             = band_label,
        area_pixel_units = round(area, 2),
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(results[seq_len(row_n)])
}

# ── 6. Analyse multiple files ─────────────────────────────────────────────────
analyze_files <- function(tif_paths, max_bands = 2) {
  all_results <- lapply(tif_paths, function(p) {
    if (!file.exists(p)) { warning("File not found: ", p); return(NULL) }
    message("\nAnalyzing: ", basename(p), "  (max bands = ", max_bands, ")")
    analyze_image(p, max_bands = max_bands)
  })
  bind_rows(all_results)
}

# ── 7. Command-line entry point ───────────────────────────────────────────────
if (!interactive()) {
  args      <- commandArgs(trailingOnly = TRUE)
  bands_arg <- grep("^--bands$", args)
  max_bands <- if (length(bands_arg) > 0) as.integer(args[bands_arg + 1]) else 2L
  tif_files <- args[!grepl("^--", args) & !seq_along(args) %in% (bands_arg + 1)]
  
  # Expand glob patterns (useful on Windows where the shell doesn't expand them)
  tif_files <- unlist(lapply(tif_files, Sys.glob))
  
  if (length(tif_files) == 0) {
    cat("Usage: Rscript analyze_wb_plots.R file1.tif [file2.tif] [--bands 2]\n")
    quit(status = 1)
  }
  
  results  <- analyze_files(tif_files, max_bands = max_bands)
  out_path <- "western_blot_areas.csv"
  write.csv(results, out_path, row.names = FALSE)
  message("\nResults saved to: ", normalizePath(out_path))
  print(results)
}