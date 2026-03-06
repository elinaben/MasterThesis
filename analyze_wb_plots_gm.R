# Load necessary libraries
if (!require("magick", quietly = TRUE)) install.packages("magick")
if (!require("pracma", quietly = TRUE)) install.packages("pracma")
library(magick)
library(pracma)

# --- 1. Image Loading and Preprocessing ---

# Read the image
image_path <- "ImagesTif/Plots_of_ELINA2026-02-1312h40m54s_IRDye800CW_.tif"
img <- image_read(image_path)

# Convert image to a grayscale matrix
# The image is grayscale, so we can extract one channel.
# Values are 0 (black) to 255 (white).
img_data <- as.integer(image_data(img, channels = "gray"))[, , 1]
img_matrix <- t(img_data) # Transpose to match R's matrix orientation (rows x cols)

# --- 2. Lane Segmentation ---

# The lanes are separated by horizontal black lines.
# We'll find these lines to define the boundaries of each lane's plot.

# Calculate the mean intensity of each row.
# Black lines will appear as deep valleys in the profile.
row_profile <- rowMeans(img_matrix)

# Invert the profile so the black lines become peaks.
inv_profile <- 255 - row_profile

# Find the peaks corresponding to the separator lines.
# We expect lines to be significantly darker than the background.
# minpeakheight = 150 ensures we pick up the black lines.
# minpeakdistance = 50 prevents detecting the same thick line multiple times.
line_peaks <- findpeaks(inv_profile, minpeakheight = 150, minpeakdistance = 50)

# Extract the y-positions of the lines.
line_y <- sort(line_peaks[, 2])

# We expect 7 lines to define the 6 plots (one above each, one below the last).
# If detection fails, fallback to splitting the image into 6 equal heights.
if (length(line_y) != 7) {
  cat("Warning: Could not detect exactly 7 separator lines. Falling back to equal splitting.\n")
  img_height <- nrow(img_matrix)
  # Exclude top 5% for potential text header if lines aren't found
  start_y <- round(img_height * 0.05)
  line_y <- round(seq(start_y, img_height, length.out = 7))
}

# --- 3. Analyze Each Lane ---

results_df <- data.frame(
  Lane = integer(),
  Band = character(),
  Position_X = numeric(),
  Peak_Height = numeric(),
  Area = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:6) {
  # Define the rows for the current lane
  y_start <- line_y[i] + 2 # Add a small buffer to avoid the line itself
  y_end <- line_y[i+1] - 2
  
  # Extract the sub-image for this lane
  lane_img <- img_matrix[y_start:y_end, ]
  
  # --- Curve Extraction ---
  # The plot is a black curve on a white background.
  # For each column (x-position), find the row (y-position) of the darkest pixel.
  # This converts the 2D plot into a 1D signal representing the curve's shape.
  raw_signal <- apply(lane_img, 2, which.min)
  
  # Create an x-axis vector
  x <- 1:length(raw_signal)
  
  # --- Signal Smoothing ---
  # Smooth the signal to reduce noise and improve peak detection.
  # Loess smoothing is effective here.
  smooth_fit <- loess(raw_signal ~ x, span = 0.15)
  signal <- predict(smooth_fit, x)
  
  # --- Peak Detection & Quantification ---
  # Find the two largest peaks in the smoothed signal.
  # These correspond to the two bands in the Western Blot.
  # npeaks = 2: we are looking for exactly 2 bands.
  # sortstr = TRUE: sort results by peak height (descending).
  peaks <- findpeaks(signal, npeaks = 2, sortstr = TRUE, minpeakdistance = 20)
  
  # Check if 2 peaks were successfully found.
  if (is.null(peaks) || nrow(peaks) < 2) {
    cat(sprintf("Warning: Could not find 2 distinct peaks in Lane %d. Skipping.\n", i))
    next
  }
  
  # Sort peaks by position (x-coordinate) to consistently label "Band 1" (left) and "Band 2" (right).
  peaks <- peaks[order(peaks[, 2]), ]
  
  for (j in 1:2) {
    # Get peak boundaries from findpeaks output
    p_start_idx <- peaks[j, 3]
    p_end_idx <- peaks[j, 4]
    
    # Extract the segment of the signal corresponding to the peak
    peak_x <- x[p_start_idx:p_end_idx]
    peak_y <- signal[p_start_idx:p_end_idx]
    
    # --- Area Calculation ---
    # Calculate the area under the curve for the peak.
    # To account for the background, we define a local baseline connecting
    # the start and end points of the peak.
    
    # Create the baseline
    baseline <- approx(x = c(peak_x[1], peak_x[length(peak_x)]),
                       y = c(peak_y[1], peak_y[length(peak_y)]),
                       xout = peak_x)$y
    
    # Calculate the area between the signal curve and the baseline using the trapezoidal rule.
    # We integrate the difference (signal - baseline).
    area <- trapz(peak_x, peak_y - baseline)
    
    # Append results to the main data frame
    results_df <- rbind(results_df, data.frame(
      Lane = i,
      Band = paste("Band", j),
      Position_X = peaks[j, 2], # x-position of the peak center
      Peak_Height = peaks[j, 1], # y-value of the peak (from baseline of the plot box)
      Area = round(area, 2) # Quantified area
    ))
  }
}

# --- 4. Output Results ---

# Define the output file name
output_csv <- "wb_analysis_results_800.csv"

# Write the results to a CSV file
write.csv(results_df, file = output_csv, row.names = FALSE)

# Print the results to the console
cat("\nAnalysis complete. Results:\n")
print(results_df)
cat(sprintf("\nResults have been saved to '%s'.\n", output_csv))