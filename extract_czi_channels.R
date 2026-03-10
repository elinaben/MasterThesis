################################################################################
# extract_czi_channels.R
#
# Extracts fluorescence channels from .czi files as colorized TIFs:
#   - Each channel saved individually in its assigned color
#   - Full composite of all channels merged
#   - DAPI + each other channel as pairwise composites
#
# Handles both standard and pyramid (multi-resolution) CZI files.
#
# Color mapping (customize COLOR_MAP below if dye names differ):
#   DAPI        -> blue
#   AF488 / 488 -> green
#   AF568 / 555 -> red
#   dodec / 594 -> magenta
#
# Output folder: <ZeissImages>/channels_<subfolder>_<YYYYMMDD_HHMMSS>/
#
# Dependencies (auto-installed on first run):
#   BiocManager::install("RBioFormats")
#   install.packages(c("tiff", "jpeg"))
################################################################################

# ── 0. Packages ───────────────────────────────────────────────────────────────

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("RBioFormats", quietly = TRUE))
  BiocManager::install("RBioFormats", ask = FALSE)
for (pkg in c("tiff", "jpeg")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(RBioFormats)
library(tiff)
library(jpeg)

# ── 1. Configuration ──────────────────────────────────────────────────────────

INPUT_DIR     <- "./ZeissImages/MF_090825/"
OUTPUT_FORMAT <- "tif"    # "tif" (16-bit lossless) or "jpg"
LOW_PCT       <- 0.1      # best-fit lower percentile clip
HIGH_PCT      <- 99.9     # best-fit upper percentile clip
Z_PROJ        <- "max"    # Z-stack collapse: "max" | "mean" | "first"

# Color map: dye short-name (lowercase, partial match) -> RGB weights c(R, G, B)
COLOR_MAP <- list(
  "dapi"  = c(0, 0, 1),   # blue
  "488"   = c(0, 1, 0),   # green
  "af488" = c(0, 1, 0),   # green
  "555"   = c(1, 0, 0),   # red
  "af555" = c(1, 0, 0),   # red
  "568"   = c(1, 0, 0),   # red
  "af568" = c(1, 0, 0),   # red
  "594"   = c(1, 0, 1),   # magenta
  "dodec" = c(1, 0, 1)    # magenta
)

# Metadata key for fluorophore short names in Zeiss CZI files
FLUOR_KEY_BASE <- paste0(
  "Experiment|AcquisitionBlock|RegionsSetup|TilesSetup|",
  "MultiTrackSetup|Track|Channel|FluorescenceDye|ShortName"
)

# ── 2. Build timestamped output directory ─────────────────────────────────────

input_clean <- gsub("[\\/]+$", "", INPUT_DIR)
folder_name <- basename(input_clean)
zeiss_root  <- dirname(input_clean)
timestamp   <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR  <- file.path(zeiss_root,
                         paste0("channels_", folder_name, "_", timestamp))

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("Output folder: ", normalizePath(OUTPUT_DIR, mustWork = FALSE), "\n")

# ── 3. Helpers ────────────────────────────────────────────────────────────────

best_fit <- function(mat) {
  lo <- quantile(mat, probs = LOW_PCT  / 100, na.rm = TRUE)
  hi <- quantile(mat, probs = HIGH_PCT / 100, na.rm = TRUE)
  if (hi == lo) return(matrix(0, nrow(mat), ncol(mat)))
  clipped <- pmin(pmax(mat, lo), hi)
  (clipped - lo) / (hi - lo)
}

get_color <- function(dye_name) {
  key <- tolower(dye_name)
  if (!is.null(COLOR_MAP[[key]])) return(COLOR_MAP[[key]])
  for (pattern in names(COLOR_MAP)) {
    if (grepl(pattern, key, fixed = TRUE)) return(COLOR_MAP[[pattern]])
  }
  message("  WARNING: No color mapping for '", dye_name, "' — using white")
  c(1, 1, 1)
}

colorize <- function(mat_01, rgb_weights) {
  dims <- dim(mat_01)
  arr  <- array(0, dim = c(dims[1], dims[2], 3))
  arr[, , 1] <- mat_01 * rgb_weights[1]
  arr[, , 2] <- mat_01 * rgb_weights[2]
  arr[, , 3] <- mat_01 * rgb_weights[3]
  arr
}

merge_colors <- function(color_array_list) {
  pmin(Reduce("+", color_array_list), 1)
}

save_image <- function(img_01, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (OUTPUT_FORMAT == "tif") {
    tiff::writeTIFF(img_01, path, bits.per.sample = 16L, compression = "none")
  } else {
    jpeg::writeJPEG(img_01, path, quality = 95)
  }
  invisible(path)
}

slice_axis <- function(arr, axis, idx) {
  index <- vector("list", length(dim(arr)))
  for (d in seq_along(index)) index[[d]] <- TRUE
  index[[axis]] <- idx
  do.call(`[`, c(list(arr), index, list(drop = TRUE)))
}

# ── 4. Collect files ──────────────────────────────────────────────────────────

czi_files <- list.files(INPUT_DIR, pattern = "\\.czi$",
                        full.names = TRUE, ignore.case = TRUE)
if (length(czi_files) == 0) stop("No .czi files found in: ", INPUT_DIR)
message(sprintf("Found %d .czi file(s).\n", length(czi_files)))

# ── 5. Main loop ──────────────────────────────────────────────────────────────

for (czi_path in czi_files) {

  base_name <- tools::file_path_sans_ext(basename(czi_path))
  out_sub   <- file.path(OUTPUT_DIR, base_name)
  dir.create(out_sub, showWarnings = FALSE, recursive = TRUE)

  message("── Processing: ", base_name)

  # Read file
  img <- tryCatch(
    RBioFormats::read.image(czi_path),
    error = function(e) { message("  ERROR reading file: ", e$message); NULL }
  )
  if (is.null(img)) next

  # ── Handle pyramid (multi-resolution) CZI ──────────────────────────────
  if (inherits(img, "AnnotatedImageList")) {
    message("  Pyramid CZI detected — using full-resolution series 1")
    img <- img[[1]]
  }

  # ── Array dimensions ────────────────────────────────────────────────────
  arr         <- as.array(img)
  actual_dims <- dim(arr)

  # ── Metadata (may be empty for pyramid files — fall back to dim()) ──────
  meta <- tryCatch(RBioFormats::coreMetadata(img), error = function(e) NULL)

  n_z <- if (!is.null(meta) && !is.null(meta$sizeZ) &&
              length(meta$sizeZ) > 0 && !is.na(meta$sizeZ)) {
    meta$sizeZ
  } else 1L

  n_channels <- if (!is.null(meta) && !is.null(meta$sizeC) &&
                    length(meta$sizeC) > 0 && !is.na(meta$sizeC)) {
    meta$sizeC
  } else {
    # For [X, Y, C] layout the channel count is the last (smallest) dim
    actual_dims[length(actual_dims)]
  }

  message(sprintf("  Array dims: [%s]  C=%d  Z=%d",
                  paste(actual_dims, collapse = " x "), n_channels, n_z))

  # ── Channel axis detection ───────────────────────────────────────────────
  c_cand <- which(actual_dims == n_channels)
  c_cand <- c_cand[c_cand > 2]
  if (length(c_cand) == 0) c_cand <- which(actual_dims == n_channels)

  if (length(c_cand) == 0 || is.na(c_cand[1])) {
    message("  ERROR: Could not detect channel axis — skipping.")
    next
  }
  c_axis <- c_cand[1]

  # ── Z axis detection ─────────────────────────────────────────────────────
  z_axis <- NULL
  if (n_z > 1) {
    z_cand <- which(actual_dims == n_z)
    z_cand <- z_cand[z_cand > 2 & z_cand != c_axis]
    if (length(z_cand) > 0) z_axis <- z_cand[1]
  }

  # ── Channel names from globalMetadata ────────────────────────────────────
  ch_names <- tryCatch({
    global <- globalMetadata(img)
    sapply(seq_len(n_channels), function(i) {
      nm <- global[[sprintf("%s #%d", FLUOR_KEY_BASE, i)]]
      if (is.null(nm) || nchar(trimws(nm)) == 0)
        sprintf("CH%02d", i)
      else
        gsub("[^[:alnum:]_.-]", "_", trimws(nm))
    })
  }, error = function(e) sprintf("CH%02d", seq_len(n_channels)))

  message("  Channels: ", paste(ch_names, collapse = ", "))

  # ── Extract, stretch, colorize ───────────────────────────────────────────
  gray_mats  <- vector("list", n_channels)
  color_arrs <- vector("list", n_channels)
  dapi_idx   <- NA

  for (i in seq_len(n_channels)) {

    ch_arr <- slice_axis(arr, c_axis, i)

    if (!is.null(z_axis) && length(dim(ch_arr)) >= 3) {
      z_now   <- if (z_axis > c_axis) z_axis - 1L else z_axis
      spatial <- setdiff(seq_len(length(dim(ch_arr))), z_now)
      ch_mat  <- switch(Z_PROJ,
        "max"   = apply(ch_arr, spatial, max),
        "mean"  = apply(ch_arr, spatial, mean),
        "first" = slice_axis(ch_arr, z_now, 1)
      )
    } else {
      ch_mat <- ch_arr
    }

    ch_mat  <- matrix(as.numeric(ch_mat), nrow = nrow(ch_mat))
    gray    <- best_fit(ch_mat)
    rgb_w   <- get_color(ch_names[i])
    col_arr <- colorize(gray, rgb_w)

    gray_mats[[i]]  <- gray
    color_arrs[[i]] <- col_arr

    if (grepl("dapi", tolower(ch_names[i]))) dapi_idx <- i

    fname <- sprintf("%s_%s.%s", base_name, ch_names[i], OUTPUT_FORMAT)
    save_image(col_arr, file.path(out_sub, fname))
    message(sprintf("  [%d/%d] %-12s  ->  %s", i, n_channels, ch_names[i], fname))
  }

  # ── Full composite ────────────────────────────────────────────────────────
  full_comp <- merge_colors(color_arrs)
  fname_comp <- sprintf("%s_composite_all.%s", base_name, OUTPUT_FORMAT)
  save_image(full_comp, file.path(out_sub, fname_comp))
  message("  Composite (all)  ->  ", fname_comp)

  # ── DAPI + each channel pairwise ─────────────────────────────────────────
  if (!is.na(dapi_idx)) {
    for (i in seq_len(n_channels)) {
      if (i == dapi_idx) next
      pair  <- merge_colors(list(color_arrs[[dapi_idx]], color_arrs[[i]]))
      fname <- sprintf("%s_DAPI+%s.%s", base_name, ch_names[i], OUTPUT_FORMAT)
      save_image(pair, file.path(out_sub, fname))
      message("  DAPI + ", ch_names[i], "  ->  ", fname)
    }
  } else {
    message("  (No DAPI channel found — skipping pairwise composites)")
  }

  message()
}

message("Done!  All images written to:\n  ",
        normalizePath(OUTPUT_DIR, mustWork = FALSE))
