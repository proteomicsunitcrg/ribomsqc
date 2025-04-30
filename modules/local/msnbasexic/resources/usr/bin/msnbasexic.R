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
json_metrics <- c("Log2_Total_Area", "dmz_ppm", "dRT_sec", "FWHM", "PPP")

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

# Create output directory if it does not exist
if (!is.null(opt$output_dir) && !dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# Read analyte data
df <- read.delim(opt$tsv_name, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Determine analytes to process
analyte_names <- if (opt$analyte_name == "all") df$short_name else c(opt$analyte_name)

# Load mzML data
ms_data <- readMSData(opt$file_name, mode = "onDisk")

# Define function to extract analyte-specific parameters
get_analyte_data <- function(df, analyte) {
  analyte_data <- subset(df, short_name == analyte)
  if (nrow(analyte_data) == 0) stop(paste("Error: Analyte not found:", analyte))
  list(mz1 = analyte_data$mz_M0, ms2_mz = analyte_data$ms2_mz, rt_target = analyte_data$rt_teoretical)
}

# Assign correct m/z based on MS level
assign_mz1 <- function(analyte_data, msLevel) {
  if (msLevel == 1) as.numeric(analyte_data$mz1) else as.numeric(analyte_data$ms2_mz)
}

# Extract chromatographic signal (XIC)
compute_xic <- function(ms_data, mz1, rt_target, rt_tol_sec, mz_tol, msLevel) {
  if (msLevel == 1) {
    chromatogram(ms_data, mz = c(mz1 - mz_tol, mz1 + mz_tol), rt = c(rt_target - rt_tol_sec, rt_target + rt_tol_sec), msLevel = 1, aggregationFun = "mean")
  } else {
    ms2_data <- filterMsLevel(ms_data, msLevel = 2)
    chromatogram(ms2_data, mz = c(mz1 - mz_tol, mz1 + mz_tol), rt = c(rt_target - rt_tol_sec, rt_target + rt_tol_sec), msLevel = 2, aggregationFun = "sum")
  }
}

# Extract max intensity m/z at observed retention time
find_mz_for_intensity_obs <- function(ms_data, mz, rt_obs, intensity_mz1) {
  spectra_near_rt <- filterRt(ms_data, rt = c(rt_obs - 1, rt_obs + 1))
  if (length(spectra_near_rt) == 0) stop("No spectra found near RT")
  all_mz_exp <- all_intensities <- c()
  for (i in seq_along(spectra_near_rt)) {
    spec <- spectra(spectra_near_rt)[[i]]
    valid_idx <- mz(spec) >= mz[1] & mz(spec) <= mz[2]
    all_mz_exp <- c(all_mz_exp, mz(spec)[valid_idx])
    all_intensities <- c(all_intensities, intensity(spec)[valid_idx])
  }
  if (length(all_mz_exp) == 0) stop("No m/z values in range")
  df <- data.frame(mz_exp = all_mz_exp, intensity = all_intensities)
  df[which.min(abs(df$intensity - intensity_mz1)), "mz_exp"]
}

# Trapezoidal area calculation
compute_xic_area <- function(rt_values, intensities) {
  valid <- which(!is.na(rt_values) & !is.na(intensities))
  if (length(valid) > 1) trapz(rt_values[valid], intensities[valid]) else NA
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
  if (length(rt_values) < 2 || all(is.na(intensities))) return(NA)
  half_max <- max(intensities, na.rm = TRUE) / 2
  idx <- which(intensities >= half_max)
  if (length(idx) < 2) return(NA)
  rt_values[max(idx)] - rt_values[min(idx)]
}

# Plot chromatogram
plot_xic <- function(rt_values, intensities, output_file) {
  if (length(rt_values) == 0 || length(intensities) == 0) return(NULL)
  xic_df <- data.frame(RT = rt_values, Intensity = intensities)
  p <- ggplot(xic_df, aes(x = RT, y = Intensity)) + geom_point(size = 1.6) + theme_minimal()
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
  analyte_data <- get_analyte_data(df, analyte)
  mz1 <- assign_mz1(analyte_data, opt$msLevel)
  mz_tol <- (opt$mz_tol_ppm * mz1) / 1e6
  rt_target <- analyte_data$rt_target
  xic <- compute_xic(ms_data, mz1, rt_target, opt$rt_tol_sec, mz_tol, opt$msLevel)
  if (length(xic) == 0 || is.null(xic[[1]])) next
  rt_values <- rtime(xic[[1]])
  int <- intensity(xic[[1]])

  if (length(rt_values) == 0 || length(int) == 0 || all(is.na(int))) {
    warning(glue("No valid XIC signal for analyte {analyte}. Skipping."))
    next
  }

  max_idx <- which.max(int)
  rt_obs <- rt_values[max_idx]
  intensity_mz1 <- int[max_idx]
  ppp_value <- calculate_ppp_custom(rt_values, int, threshold_ratio = 0.01)
  fwhm_value <- calculate_fwhm(rt_values, int)
  area_mz1 <- compute_xic_area(rt_values, int)
  mz_exp_result <- find_mz_for_intensity_obs(ms_data, c(mz1 - mz_tol, mz1 + mz_tol), rt_obs, intensity_mz1)
  sample_base <- tools::file_path_sans_ext(basename(opt$file_name))
  sample_name <- sub("\\.raw$", "", sample_base)
  result_row <- function_create_output_df(sample_name, analyte, rt_target, rt_obs, mz1, mz_exp_result, intensity_mz1, area_mz1, fwhm_value, ppp_value)
  all_results[[analyte]] <- result_row
  for (metric in json_metrics) {
    if (metric %in% colnames(result_row)) {
      update_metric_json(metric, analyte, sample_name, result_row[[metric]])
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