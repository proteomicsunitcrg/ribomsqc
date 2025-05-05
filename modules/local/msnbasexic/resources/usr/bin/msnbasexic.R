#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(MSnbase))
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(glue))

# Define metrics to be exported to JSON
json_metrics <- c("Log2_Total_Area", "dmz_ppm", "Observed_RT_sec", "dRT_sec", "FWHM", "PPP")

# PPP constant:  
PPP_THRESHOLD_RATIO <- 0.05

# Define command-line options
option_list <- list(
  make_option(c("--file_name"), type = "character", help = "Input mzML file path"),
  make_option(c("--tsv_name"), type = "character", help = "Input TSV file with analyte info"),
  make_option(c("--output_dir"), type = "character", help = "Directory to write output files"),
  make_option(c("--analyte_name"), type = "character", help = "Name of analyte to process"),
  make_option(c("--rt_tol_sec"), type = "numeric", help = "Retention time tolerance in seconds"),
  make_option(c("--mz_tol_ppm"), type = "numeric", help = "m/z tolerance in ppm"),
  make_option(c("--msLevel"), type = "numeric", help = "MS level to use"),
  make_option(c("--plot_xic_ms1"), type = "logical", default = FALSE, help = "Plot XIC for MS1"),
  make_option(c("--plot_xic_ms2"), type = "logical", default = FALSE, help = "Plot XIC for MS2"),
  make_option(c("--plot_output_path"), type = "character", help = "Path to save XIC plot"),
  make_option(c("--overwrite_tsv"), type = "logical", default = FALSE, help = "Overwrite existing TSV")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))
opt$plot_xic_ms1 <- as.logical(opt$plot_xic_ms1)
opt$plot_xic_ms2 <- as.logical(opt$plot_xic_ms2)
opt$overwrite_tsv <- as.logical(opt$overwrite_tsv)

# Create output directory if it does not exist
if (!is.null(opt$output_dir) && !dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
  message(glue("Created output directory: {opt$output_dir}"))
}

# Read analyte data
df <- read.delim(opt$tsv_name, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
message(glue("Read {nrow(df)} analytes from TSV: {opt$tsv_name}"))

# Determine analytes to process
analyte_names <- if (opt$analyte_name == "all") df$short_name else c(opt$analyte_name)
message(glue("Analytes to process: {paste(analyte_names, collapse = ', ')}"))

# Load mzML data
message(glue("Loading mzML file: {opt$file_name}"))
ms_data <- readMSData(opt$file_name, mode = "onDisk")
message("MS data loaded.")

# Define function to extract analyte-specific parameters
get_analyte_data <- function(df, analyte) {
  analyte_data <- subset(df, short_name == analyte)
  if (nrow(analyte_data) == 0) stop(paste("Error: Analyte not found:", analyte))
  list(mz1 = analyte_data$mz_M0, ms2_mz = analyte_data$ms2_mz, rt_target = analyte_data$rt_teoretical)
}

# Assign m/z value based on MS level
assign_mz1 <- function(analyte_data, msLevel) {
  if (msLevel == 1) {
    as.numeric(analyte_data$mz1)
  } else {
    as.numeric(analyte_data$mz1)  # Use precursor m/z for MS2 as well
  }
}

# Compute XIC for given m/z and retention time
compute_xic <- function(ms_data, mz1, rt_target, rt_tol_sec, mz_tol, msLevel, ppm_tol = 20) {
  rt_window <- c(rt_target - rt_tol_sec, rt_target + rt_tol_sec)
  mz_window <- c(mz1 - mz_tol, mz1 + mz_tol)

  if (msLevel == 1) {
    # MS1: extract chromatogram directly using mz and RT window
    chromatogram(ms_data, mz = mz_window, rt = rt_window, msLevel = 1, aggregationFun = "mean")
  } else {
    # MS2: filter manually by RT and precursor m/z
    ms2_data <- ms_data |>
      filterMsLevel(2) |>
      filterRt(rt_window)

    # Extract precursor m/z for each spectrum
    prec_mz <- precursorMz(ms2_data)

    # Keep only spectra with non-NA precursor m/z within target range
    valid_idx <- which(!is.na(prec_mz) & prec_mz >= mz_window[1] & prec_mz <= mz_window[2])

    # Logging: how many MS2 spectra matched
    message(glue("Found {length(valid_idx)} MS2 spectra with precursor m/z within ±{round(ppm_tol, 1)} ppm of {mz1}"))

    if (length(valid_idx) == 0) {
      warning("No MS2 spectra matched precursor m/z within the specified ppm.")
      return(list(NULL))
    }

    # Keep only matching MS2 spectra
    ms2_data <- ms2_data[valid_idx]

    # Build chromatogram manually by summing intensities per spectrum
    rt_vals <- rtime(ms2_data)
    int_vals <- sapply(spectra(ms2_data), function(s) sum(intensity(s), na.rm = TRUE))

    chrom <- new("Chromatogram", rtime = rt_vals, intensity = int_vals)
    list(chrom)
  }
}

# Extract max intensity m/z at observed retention time
find_mz_for_intensity_obs <- function(ms_data, mz, rt_obs, intensity_mz1, msLevel) {
  data_filtered <- ms_data |>
    filterMsLevel(msLevel) |>
    filterRt(rt = c(rt_obs - 1, rt_obs + 1))

  if (length(data_filtered) == 0) stop("No spectra found near RT")

  all_mz_exp <- all_intensities <- c()

  for (spec in spectra(data_filtered)) {
    valid_idx <- mz(spec) >= mz[1] & mz(spec) <= mz[2]
    all_mz_exp <- c(all_mz_exp, mz(spec)[valid_idx])
    all_intensities <- c(all_intensities, intensity(spec)[valid_idx])
  }

  if (length(all_mz_exp) == 0) stop("No m/z values in range")

  df <- data.frame(mz_exp = all_mz_exp, intensity = all_intensities)
  df[which.min(abs(df$intensity - intensity_mz1)), "mz_exp"]
}

# Format numeric values for consistency
format_custom <- function(x) {
  if (is.na(x)) return("NA")
  if (abs(x) >= 1e4) format(x, scientific = TRUE, digits = 3) else format(x, scientific = FALSE, digits = 3)
}

# Construct results row for one analyte
function_create_output_df <- function(sample_name, analyte, rt_target, rt_obs, mz1, mz_exp_result, intensity_mz1, area_mz1, fwhm, ppp) {
  data.frame(
    Sample = sample_name,
    short_name = analyte,
    Expected_RT_sec = rt_target,
    Observed_RT_sec = rt_obs,
    dRT_sec = rt_obs - rt_target,
    Expected_MZ = mz1,
    Observed_MZ = mz_exp_result,
    dmz_ppm = ((mz_exp_result - mz1) / mz1) * 1e6,
    dmz_Da = mz_exp_result - mz1,
    Max_Intensity = format_custom(intensity_mz1),
    Area = format_custom(area_mz1),
    Log2_Total_Area = format_custom(ifelse(!is.na(area_mz1) & area_mz1 > 0, log2(area_mz1), NA)),
    FWHM = format_custom(fwhm),
    PPP = format_custom(ppp)
  )
}

# Additional metric calculations
calculate_ppp_custom <- function(rt_values, intensities, threshold_ratio = 0.5) {
  max_int <- max(intensities, na.rm = TRUE)
  sum(!is.na(intensities) & intensities >= max_int * threshold_ratio)
}

calculate_fwhm <- function(rt_values, intensities) {
  if (length(rt_values) < 2 || all(is.na(intensities))) return(list(fwhm = NA, left_rt = NA, right_rt = NA))

  half_max <- max(intensities, na.rm = TRUE) / 2

  ord <- order(rt_values)
  rt_values <- rt_values[ord]
  intensities <- intensities[ord]

  peak_idx <- which.max(intensities)

  left_rt <- NA
  for (i in seq(peak_idx, 2, by = -1)) {
    if (intensities[i - 1] < half_max && intensities[i] >= half_max) {
      x1 <- rt_values[i - 1]
      x2 <- rt_values[i]
      y1 <- intensities[i - 1]
      y2 <- intensities[i]
      left_rt <- x1 + (half_max - y1) * (x2 - x1) / (y2 - y1)
      break
    }
  }

  right_rt <- NA
  for (i in seq(peak_idx, length(intensities) - 1)) {
    if (intensities[i + 1] < half_max && intensities[i] >= half_max) {
      x1 <- rt_values[i]
      x2 <- rt_values[i + 1]
      y1 <- intensities[i]
      y2 <- intensities[i + 1]
      right_rt <- x1 + (half_max - y1) * (x2 - x1) / (y2 - y1)
      break
    }
  }

  if (is.na(left_rt) || is.na(right_rt)) {
    return(list(fwhm = NA, left_rt = left_rt, right_rt = right_rt))
  }

  return(list(fwhm = right_rt - left_rt, left_rt = left_rt, right_rt = right_rt))
}

# Trapezoidal area calculation
compute_xic_area <- function(rt_values, intensities) {
  valid <- which(!is.na(rt_values) & !is.na(intensities))
  if (length(valid) > 1) {
    trapz(rt_values[valid], intensities[valid])
  } else {
    NA
  }
}

plot_xic <- function(rt_values, intensities, output_file, analyte_name = NULL, ms_level = NULL, mz_tol = NULL, rt_tol_sec = NULL, sample_name = NULL, mz_tol_ppm = NULL) {
  if (length(rt_values) == 0 || length(intensities) == 0) return(NULL)
  
  xic_df <- data.frame(RT = rt_values, Intensity = intensities)
  xic_df <- xic_df[complete.cases(xic_df), ]
  xic_df <- xic_df[order(xic_df$RT), ]

  max_idx <- which.max(xic_df$Intensity)
  max_rt <- xic_df$RT[max_idx]
  max_intensity <- xic_df$Intensity[max_idx]

  half_max <- max_intensity / 2
  ppp <- calculate_ppp_custom(xic_df$RT, xic_df$Intensity, threshold_ratio = PPP_THRESHOLD_RATIO)   

  # Get FWHM and edges
  fwhm_data <- calculate_fwhm(xic_df$RT, xic_df$Intensity)
  fwhm <- fwhm_data$fwhm
  left_rt <- fwhm_data$left_rt
  right_rt <- fwhm_data$right_rt

  # Create FWHM line with interpolated points
  fwhm_line <- NULL
  if (!is.na(left_rt) && !is.na(right_rt)) {
    fwhm_line <- data.frame(
      RT = seq(left_rt, right_rt, length.out = 100),
      Intensity = rep(half_max, 100)
    )
  }

  # PPP threshold
  ppp_threshold <- max_intensity * PPP_THRESHOLD_RATIO

  # Titles
  plot_title <- glue("XIC for {analyte_name} (MS{ms_level})")
  plot_subtitle <- glue("Sample: {sample_name} | mz tol: ±{mz_tol_ppm} ppm | rt tol: ±{round(rt_tol_sec, 1)}s")

  p <- ggplot(xic_df, aes(x = RT, y = Intensity)) +
    geom_line(color = "blue", linewidth = 0.8) +
    geom_point(size = 1.6) +
    geom_vline(xintercept = max_rt, linetype = "dashed", color = "red") +
    annotate("text", x = max_rt, y = max_intensity, label = "Max", vjust = -1, hjust = 1, size = 3.5) +
    annotate("text", x = Inf, y = Inf, label = glue("FWHM: {round(fwhm, 2)} sec\nPPP: {ppp}"), 
             hjust = 1.1, vjust = 1.3, size = 4, color = "black") +
    # Improved FWHM visual
    {if (!is.null(fwhm_line)) geom_line(data = fwhm_line, aes(x = RT, y = Intensity),
                                        inherit.aes = FALSE, color = "darkgreen",
                                        linetype = "dotted", linewidth = 1)} +
    # PPP threshold line
    geom_hline(yintercept = ppp_threshold, color = "gray40", linetype = "dotdash") +
    annotate("text", x = -Inf, y = ppp_threshold, label = "PPP threshold (1%)", 
             hjust = -0.1, vjust = -0.5, size = 3.5, color = "gray30") +
    ggtitle(plot_title, subtitle = plot_subtitle) +
    xlab("Retention Time (sec)") + ylab("Intensity") +
    theme_minimal()

  ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300, bg = "white")
}

# Write metric to JSON file
update_metric_json <- function(metric_name, analyte, sample_name, value) {
  json_file <- file.path(getwd(), paste0(metric_name, "_", sample_name, "_mqc.json"))
  if (file.exists(json_file)) {
    json_data <- fromJSON(json_file)
  } else {
    json_data <- list(
      id = metric_name,
      section_name = metric_name,
      description = paste(metric_name, "values across samples"),
      plot_type = "linegraph",
      pconfig = list(
        id = paste0(metric_name, "_plot"),
        title = metric_name,
        xlab = "Sample",
        ylab = paste(metric_name, "(unit)"),
        xlab_format = "category",
        showlegend = TRUE
      ),
      data = list()
    )
  }
  json_data$data[[analyte]][[sample_name]] <- value
  write_json(json_data, path = json_file, auto_unbox = TRUE, pretty = TRUE)
}

# Main processing loop
all_results <- list()
for (analyte in analyte_names) {
  message(glue("Processing analyte: {analyte}"))
  analyte_data <- get_analyte_data(df, analyte)
  mz1 <- assign_mz1(analyte_data, opt$msLevel)
  mz_tol <- (opt$mz_tol_ppm * mz1) / 1e6
  rt_target <- analyte_data$rt_target
  message(glue("  mz1 = {mz1}, rt_target = {rt_target}, mz_tol = {mz_tol}, msLevel = {opt$msLevel}"))
  xic <- compute_xic(ms_data, mz1, rt_target, opt$rt_tol_sec, mz_tol, opt$msLevel)
  if (length(xic) == 0 || is.null(xic[[1]])) {
    message(glue("  No XIC signal found for analyte {analyte}. Skipping."))
    next
  }
  rt_values <- rtime(xic[[1]])
  int <- intensity(xic[[1]])

  if (length(rt_values) == 0 || length(int) == 0 || all(is.na(int))) {
    warning(glue("No valid XIC signal for analyte {analyte}. Skipping."))
    next
  }

  max_idx <- which.max(int)
  rt_obs <- rt_values[max_idx]
  intensity_mz1 <- int[max_idx]
  ppp_value <- calculate_ppp_custom(rt_values, int, threshold_ratio = PPP_THRESHOLD_RATIO)
  fwhm_data <- calculate_fwhm(rt_values, int)
  fwhm_value <- fwhm_data$fwhm
  area_mz1 <- compute_xic_area(rt_values, int)
  
  mz_exp_result <- tryCatch({
    find_mz_for_intensity_obs(ms_data, c(mz1 - mz_tol, mz1 + mz_tol), rt_obs, intensity_mz1, opt$msLevel)
  }, error = function(e) {
    warning(glue("Could not determine observed m/z for {analyte}: {e$message}"))
    NA
  })

  sample_base <- tools::file_path_sans_ext(basename(opt$file_name))
  sample_name <- sub("\\.raw$", "", sample_base)
  result_row <- function_create_output_df(sample_name, analyte, rt_target, rt_obs, mz1, mz_exp_result, intensity_mz1, area_mz1, fwhm_value, ppp_value)
  all_results[[analyte]] <- result_row
  message(glue("  Calculated metrics for {analyte}:"))
  message(glue("    Observed RT: {rt_obs}"))
  message(glue("    Observed MZ: {mz_exp_result}"))
  message(glue("    Max Intensity: {intensity_mz1}"))
  message(glue("    Area: {area_mz1}"))
  message(glue("    Log2 Area: {ifelse(!is.na(area_mz1) & area_mz1 > 0, log2(area_mz1), NA)}"))
  message(glue("    FWHM: {fwhm_value}"))
  message(glue("    PPP: {ppp_value}"))

  if (!is.null(opt$plot_output_path)) {
    should_plot_ms1 <- opt$plot_xic_ms1 && opt$msLevel == 1
    should_plot_ms2 <- opt$plot_xic_ms2 && opt$msLevel == 2
    if (should_plot_ms1 || should_plot_ms2) {
      dir.create(dirname(opt$plot_output_path), recursive = TRUE, showWarnings = FALSE)
      message(glue("  Saving XIC plot to: {opt$plot_output_path}"))
      plot_xic(rt_values, int, opt$plot_output_path,
               analyte_name = analyte,
               ms_level = opt$msLevel,
               mz_tol = mz_tol,
               rt_tol_sec = opt$rt_tol_sec,
               sample_name = sample_name,
               mz_tol_ppm = opt$mz_tol_ppm)
    }
  }

  for (metric in json_metrics) {
    if (metric %in% colnames(result_row)) {
      update_metric_json(metric, analyte, sample_name, result_row[[metric]])
      message(glue("    Updated JSON metric: {metric}"))
    }
  }
}

# Combine all results into final dataframe
output_df <- do.call(rbind, all_results)

# Save version information
version_info <- list(
  msnbasexic_script = "msnbasexic.R",
  generated_at = as.character(Sys.time()),
  R_version = paste(R.version$major, R.version$minor, sep = "."),
  MSnbase_version = as.character(packageVersion("MSnbase")),
  ggplot2_version = as.character(packageVersion("ggplot2")),
  pracma_version = as.character(packageVersion("pracma"))
)
lines <- unlist(lapply(names(version_info), function(k) paste0(k, ": ", version_info[[k]])))
writeLines(lines, "versions.yml")

print("versions.yml generated.")
