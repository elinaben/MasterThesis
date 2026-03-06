# analyze_wb_plots.R  —  version 7
# ─────────────────────────────────────────────────────────────────────────────
# Requirements:  install.packages(c("tiff", "dplyr"))
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tiff)
  library(dplyr)
})

message("analyze_wb_plots.R  version 7  loaded")


# ── Peak finder: rolling-window maximum ───────────────────────────────────────
# For each position i, it is a peak if y[i] is the maximum in a window of
# radius `half_w` AND it is at least as large as immediate left/right neighbors.
# This is the simplest possible formulation, with no diff(), sign(), or ==.
# Multiple adjacent positions with the same maximum value are collapsed to the
# leftmost one, then the midpoint of the run is taken.
#
find_peaks_r <- function(y, half_w, min_prominence = 0) {
  n <- length(y)
  if (n < 3L) return(integer(0))

  is_peak <- logical(n)
  for (i in seq_len(n)) {
    lo   <- max(1L, i - half_w)
    hi   <- min(n,  i + half_w)
    win  <- y[lo:hi]
    # i is a local maximum in the window AND not lower than both immediate neighbors
    if (y[i] >= max(win) &&
        (i == 1L || y[i] >= y[i - 1L]) &&
        (i == n  || y[i] >= y[i + 1L])) {
      is_peak[i] <- TRUE
    }
  }

  if (!any(is_peak)) return(integer(0))

  # Collapse adjacent peak positions to midpoints of each run
  peak_pos <- which(is_peak)
  if (length(peak_pos) == 0L) return(integer(0))

  # Group consecutive indices, return midpoint of each group
  groups  <- list()
  grp_cur <- peak_pos[1L]
  for (k in seq_along(peak_pos)[-1L]) {
    if (peak_pos[k] == tail(grp_cur, 1L) + 1L) {
      grp_cur <- c(grp_cur, peak_pos[k])
    } else {
      groups  <- c(groups, list(grp_cur))
      grp_cur <- peak_pos[k]
    }
  }
  groups <- c(groups, list(grp_cur))

  peaks <- sapply(groups, function(g) (g[1L] + g[length(g)]) %/% 2L)

  if (length(peaks) == 0L) return(integer(0))

  # Prominence filter
  keep <- logical(length(peaks))
  for (ii in seq_along(peaks)) {
    p  <- peaks[ii]
    hl <- peaks[peaks  < p & y[peaks] > y[p]]
    hr <- peaks[peaks  > p & y[peaks] > y[p]]
    lb <- if (length(hl) > 0L) max(hl) else 1L
    rb <- if (length(hr) > 0L) min(hr) else n
    prom <- y[p] - max(min(y[lb:p]), min(y[p:rb]))
    keep[ii] <- prom >= min_prominence
  }

  peaks[keep]
}


# ── Inline test (runs on load, prints result) ─────────────────────────────────
local({
  ty <- c(1, 3, 5, 7, 9, 9, 9, 7, 5, 3, 1, 5, 8, 12, 15, 12, 8, 5, 1)
  hw <- max(floor(19L / 6L), 5L)   # = 5
  tp <- find_peaks_r(ty, half_w = hw, min_prominence = 2)
  message(sprintf("  inline test: half_w=%d  peaks=[%s]  (expect [6,15])",
                  hw, paste(tp, collapse = ",")))
})


# ── Subplot regions ────────────────────────────────────────────────────────────
find_subplot_regions <- function(arr, min_height = 20L) {
  white_rows <- rowMeans(arr) > (250 / 255)
  regions <- list(); in_plot <- FALSE; start <- NULL
  for (i in seq_along(white_rows)) {
    if (!white_rows[i] && !in_plot) {
      in_plot <- TRUE; start <- i
    } else if (white_rows[i] && in_plot) {
      in_plot <- FALSE
      if ((i - 1L - start) > min_height)
        regions <- c(regions, list(c(start, i - 1L)))
    }
  }
  if (in_plot && (nrow(arr) - start) > min_height)
    regions <- c(regions, list(c(start, nrow(arr))))
  regions
}


# ── Trace curve ────────────────────────────────────────────────────────────────
trace_curve <- function(subplot, top_skip = 5L, bottom_skip = 3L,
                        side_skip = 3L, dark_threshold = 10 / 255) {
  h        <- nrow(subplot)
  interior <- subplot[(top_skip + 1L):(h - bottom_skip),
                      (side_skip + 1L):(ncol(subplot) - side_skip)]
  n_cols  <- ncol(interior)
  curve_y <- rep(NA_real_, n_cols)
  for (col in seq_len(n_cols)) {
    dark_rows <- which(interior[, col] < dark_threshold)
    if (length(dark_rows) > 0L)
      curve_y[col] <- mean(dark_rows)
  }
  valid <- which(!is.na(curve_y))
  if (length(valid) == 0L) return(NULL)
  list(col_idx = valid, y_vals = curve_y[valid])
}


# ── Area ───────────────────────────────────────────────────────────────────────
area_under_segment <- function(y_vals) {
  n        <- length(y_vals)
  baseline <- seq(y_vals[1L], y_vals[n], length.out = n)
  signal   <- pmax(y_vals - baseline, 0)
  sum((head(signal, -1L) + tail(signal, -1L)) / 2)
}


# ── Band detection ─────────────────────────────────────────────────────────────
detect_bands <- function(y_vals, max_bands = 2L, lane_idx = NA) {
  n      <- length(y_vals)
  half_w <- max(floor(n / 6L), 5L)

  message(sprintf(
    "    [debug] lane %s: n=%d  half_w=%d  y[1..5]=[%s]  y_max=%.1f",
    lane_idx, n, half_w,
    paste(round(y_vals[1:min(5L,n)], 1), collapse = ","),
    max(y_vals)
  ))

  peaks <- find_peaks_r(y_vals, half_w = half_w, min_prominence = 5)

  message(sprintf(
    "    [debug] lane %s: peaks_found=%d  cols=[%s]  y=[%s]",
    lane_idx, length(peaks),
    paste(peaks, collapse = ","),
    paste(round(y_vals[peaks], 1), collapse = ",")
  ))

  if (length(peaks) < 2L) return(list(c(1L, n)))

  top_peaks <- sort(
    peaks[order(y_vals[peaks], decreasing = TRUE)[seq_len(min(max_bands, length(peaks)))]]
  )

  splits <- c(1L)
  for (i in seq_len(length(top_peaks) - 1L)) {
    seg_y  <- y_vals[top_peaks[i]:top_peaks[i + 1L]]
    valley <- top_peaks[i] + which.min(seg_y) - 1L
    splits <- c(splits, valley)
  }
  splits <- c(splits, n)

  lapply(seq_len(length(splits) - 1L), function(i) c(splits[i], splits[i + 1L]))
}


# ── Analyse one TIFF ───────────────────────────────────────────────────────────
analyze_image <- function(tif_path, max_bands = 2L) {
  arr <- readTIFF(tif_path, as.is = FALSE)
  if (length(dim(arr)) == 3L)
    arr <- arr[,,1] * 0.299 + arr[,,2] * 0.587 + arr[,,3] * 0.114

  regions <- find_subplot_regions(arr)
  message(sprintf("  %s: %d lanes detected", basename(tif_path), length(regions)))

  rows <- vector("list", length(regions) * max_bands); row_n <- 0L

  for (lane_idx in seq_along(regions)) {
    r       <- regions[[lane_idx]]
    subplot <- arr[r[1]:r[2], ]
    curve   <- trace_curve(subplot)
    if (is.null(curve)) { message(sprintf("  [WARN] lane %d: no curve", lane_idx)); next }

    bands <- detect_bands(curve$y_vals, max_bands = max_bands, lane_idx = lane_idx)

    for (band_idx in seq_along(bands)) {
      seg        <- bands[[band_idx]]
      segment    <- curve$y_vals[seg[1]:seg[2]]
      area       <- area_under_segment(segment)
      band_label <- if (length(bands) > 1L) paste0("band_", band_idx) else "single"
      message(sprintf("  Lane %d | %s: area = %.2f px\u00b2", lane_idx, band_label, area))
      row_n <- row_n + 1L
      rows[[row_n]] <- data.frame(image = basename(tif_path), lane = lane_idx,
                                  band = band_label, area_pixel_units = round(area, 2L),
                                  stringsAsFactors = FALSE)
    }
  }
  bind_rows(rows[seq_len(row_n)])
}


# ── Analyse multiple files ─────────────────────────────────────────────────────
analyze_files <- function(tif_paths, max_bands = 2L) {
  all_results <- lapply(tif_paths, function(p) {
    if (!file.exists(p)) { warning("File not found: ", p); return(NULL) }
    message("\nAnalyzing: ", basename(p), "  (max bands = ", max_bands, ")")
    results  <- analyze_image(p, max_bands = max_bands)
    out_path <- file.path(dirname(p),
                          paste0(tools::file_path_sans_ext(basename(p)), "_quantification.csv"))
    write.csv(results, out_path, row.names = FALSE)
    message("  Saved: ", normalizePath(out_path, mustWork = FALSE))
    results
  })
  bind_rows(all_results)
}


# ── CLI ────────────────────────────────────────────────────────────────────────
if (!interactive()) {
  args      <- commandArgs(trailingOnly = TRUE)
  bands_arg <- which(args == "--bands")
  max_bands <- if (length(bands_arg) > 0L) as.integer(args[bands_arg + 1L]) else 2L
  drop_idx  <- c(bands_arg, bands_arg + 1L)
  tif_files <- if (length(drop_idx) > 0L) args[-drop_idx] else args
  tif_files <- unlist(lapply(tif_files, Sys.glob))
  if (length(tif_files) == 0L) {
    cat("Usage: Rscript analyze_wb_plots.R file1.tif [file2.tif] [--bands 2]\n")
    quit(status = 1L)
  }
  results <- analyze_files(tif_files, max_bands = max_bands)
  print(results)
}
