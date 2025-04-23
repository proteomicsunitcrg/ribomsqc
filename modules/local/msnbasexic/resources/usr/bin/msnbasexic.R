#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(MSnbase))
suppressPackageStartupMessages(library(pracma))  
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--file_name"), type = "character", help = "Input mzML file path"),
  make_option(c("--tsv_name"), type = "character", help = "Input TSV file with analyte info"),
  make_option(c("--output_dir"), type = "character", help = "Directory for output TSV file"),
  make_option(c("--analyte_name"), type = "character", help = "Name of analyte to process"),
  make_option(c("--rt_tol_sec"), type = "numeric", help = "Retention time tolerance in seconds"),
  make_option(c("--mz_tol_ppm"), type = "numeric", help = "m/z tolerance in ppm"),
  make_option(c("--msLevel"), type = "numeric", help = "MS level to use"),
  make_option(c("--plot_xic_ms1"), type = "logical", default = FALSE, help = "Plot XIC for MS1"),
  make_option(c("--plot_xic_ms2"), type = "logical", default = FALSE, help = "Plot XIC for MS2"),
  make_option(c("--plot_output_path"), type = "character", help = "Path to save XIC plot"),
  make_option(c("--overwrite_tsv"), type = "logical", default = FALSE, help = "Overwrite existing TSV")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Extracts m/z values and retention times from the TSV file.
get_analyte_data <- function(df, analyte) {
  analyte_data <- subset(df, short_name == analyte)

  if (nrow(analyte_data) == 0) {
    stop(paste("Error: Analyte not found:", analyte))
  }

  return(list(
    mz1 = analyte_data$mz_M0,
    ms2_mz = analyte_data$ms2_mz,
    rt_target = analyte_data$rt_teoretical
  ))
}

# Determines which m/z to use for XIC extraction depending on MS1 or MS2
assign_mz1 <- function(analyte_data, msLevel) {
  if (msLevel == 1) {
    print(paste("MS level 1 with mass", analyte_data$mz1))
    return(as.numeric(analyte_data$mz1))
  } else {
    print(paste("MS level 2 with transition mass", analyte_data$ms2_mz))
    return(as.numeric(analyte_data$ms2_mz))
  }
}

#  Extracts the chromatogram for a given m/z and retention time window.
compute_xic <- function(ms_data, mz1, rt_target, rt_tol_sec, mz_tol, msLevel) {
  print("Computing XIC...")

  if (msLevel == 1) {
    xic <- chromatogram(ms_data, 
                        mz = c(mz1 - mz_tol, mz1 + mz_tol), 
                        rt = c(rt_target - rt_tol_sec, rt_target + rt_tol_sec), 
                        msLevel = 1, 
                        aggregationFun = "mean")
  } else {
    ms2_data <- filterMsLevel(ms_data, msLevel = 2)
    print(paste("MS2 spectra after RT filtering:", length(ms2_data)))

    xic <- chromatogram(ms2_data, 
                        mz = c(mz1 - mz_tol, mz1 + mz_tol), 
                        rt = c(rt_target - rt_tol_sec, rt_target + rt_tol_sec), 
                        msLevel = 2, 
                        aggregationFun = "sum")
  }

  return(xic)
}

#  Finds the closest m/z value at the retention time of the highest intensity.
find_mz_for_intensity_obs <- function(ms_data, mz, rt_obs, intensity_mz1) {
  spectra_near_rt <- filterRt(ms_data, rt = c(rt_obs - 1, rt_obs + 1))  # ±1s
  
  if (length(spectra_near_rt) == 0) {
    stop("Error: No spectra found near the observed RT.")
  }

  all_mz_exp <- c()
  all_intensities <- c()

  for (i in seq_along(spectra_near_rt)) {
    spec <- spectra(spectra_near_rt)[[i]]
    mz_spec <- mz(spec)
    intensity_spec <- intensity(spec)

    valid_idx <- mz_spec >= mz[1] & mz_spec <= mz[2]

    if (sum(valid_idx) > 0) {
      all_mz_exp <- c(all_mz_exp, mz_spec[valid_idx])
      all_intensities <- c(all_intensities, intensity_spec[valid_idx])
    }
  }

  if (length(all_mz_exp) == 0) {
    stop("Error: No m/z values found within the specified range.")
  }

  df_mz_exp <- data.frame(mz_exp = all_mz_exp, intensity = all_intensities)
  closest_match <- df_mz_exp[which.min(abs(df_mz_exp$intensity - intensity_mz1)), ]

  return(closest_match$mz_exp)
}

# Computes the trapezoidal area of the XIC.
compute_xic_area <- function(rt_values, intensities) {
  valid_idx <- which(!is.na(rt_values) & !is.na(intensities))

  if (length(valid_idx) > 1) {
    area <- trapz(rt_values[valid_idx], intensities[valid_idx])
  } else {
    area <- NA
  }

  return(area)
}

# Ensures numerical values are formatted correctly.
format_custom <- function(x) {
  if (is.na(x)) {
    return("NA")
  } else if (abs(x) >= 1e4) {
    return(format(x, scientific = TRUE, digits = 3))
  } else {
    return(format(x, scientific = FALSE, digits = 3))
  }
}

#  Organizes extracted results into a structured table. 
function_create_output_df <- function(sample_name, analyte, rt_target, rt_obs, mz1, mz_exp_result, intensity_mz1, area_mz1, fwhm, ppp){

  output_df <- data.frame(
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

  return(output_df)
}

# Function to write results to a TSV file with MultiQC header
write_results <- function(output_df, file_name) {
  
  # Generate output filename
  sample_name <- tools::file_path_sans_ext(basename(file_name))
  output_file <- paste0(sample_name, "_mqc.tsv")
  
  # Define MultiQC header lines
  header_lines <- c(
    "# id: 'qc_table'",
    "# section_name: 'List of detected nucleosides'",
    "# plot_type: 'table'"
  )
  
  # Open a connection to the output file
  con <- file(output_file, open = "wt")
  
  # Write header lines
  writeLines(header_lines, con)
  
  # Write a blank line to separate header from table (important!)
  writeLines("", con)
  
  # Write the data frame as TSV
  write.table(output_df, con, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Close the file connection
  close(con)
  
  print(paste("Results saved in:", output_file))
}

calculate_ppp_custom <- function(rt_values, intensities, threshold_ratio = 0.5) {
  if (length(intensities) == 0 || all(is.na(intensities))) return(NA)

  max_int <- max(intensities, na.rm = TRUE)
  threshold <- max_int * threshold_ratio

  valid_idx <- which(!is.na(intensities) & intensities >= threshold)
  return(length(valid_idx))
}

calculate_fwhm <- function(rt_values, intensities) {
  if (length(rt_values) < 2 || length(intensities) < 2) return(NA)
  if (all(is.na(intensities)) || max(intensities, na.rm = TRUE) == 0) return(NA)

  max_int <- max(intensities, na.rm = TRUE)
  half_max <- max_int / 2

  valid_idx <- which(!is.na(intensities) & intensities >= half_max)

  if (length(valid_idx) < 2) return(NA)

  rt_left <- rt_values[min(valid_idx)]
  rt_right <- rt_values[max(valid_idx)]

  return(rt_right - rt_left)
}

plot_xic <- function(rt_values, intensities, output_file) {
  if (length(rt_values) == 0 || length(intensities) == 0) {
    warning("Warning: No data to plot.")
    return(NULL)
  }

  show_line <- FALSE  # ← Internal control

  plot_dir <- dirname(output_file)
  if (plot_dir == ".") {
    plot_dir <- getwd()
  }
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  xic_df <- data.frame(
    RT = rt_values,
    Intensity = intensities,
    Signal = ifelse(!is.na(intensities) & intensities > 0, "PPP", "Zero/NA")
  ) |> na.omit()

  xic_plot <- ggplot(xic_df, aes(x = RT, y = Intensity)) +
    scale_color_manual(values = c("PPP" = "blue", "Zero/NA" = "gray80")) +
    labs(title = "Extracted Ion Chromatogram (XIC)",
         subtitle = paste("PPP =", sum(xic_df$Signal == "PPP")),
         x = "Retention Time (s)",
         y = "Intensity") +
    theme_minimal()

  if (show_line) {
    xic_plot <- xic_plot + geom_line(linewidth = 0.5, color = "black", alpha = 0.6)
  }

  xic_plot <- xic_plot + geom_point(aes(color = Signal), size = 1.6, alpha = 0.9)

  ggsave(output_file, plot = xic_plot, width = 8, height = 6, dpi = 300, bg = "white")
  cat("XIC plot saved:", output_file, "\n")
}


setwd(getwd())  # Ensure working directory is correctly set
plot_xic_ms1 <- opt$plot_xic_ms1
plot_xic_ms2 <- opt$plot_xic_ms2
plot_output_path <- opt$plot_output_path

if (is.null(opt$file_name) || is.null(opt$tsv_name) || is.null(opt$analyte_name) ||
    is.null(opt$rt_tol_sec) || is.null(opt$mz_tol_ppm) || is.null(opt$msLevel) ||
    is.null(opt$output_dir)) {
  stop("Error: Missing required parameters. Ensure all required inputs are provided.")
}



print(paste("Analyte:", opt$analyte_name))
df <- read.delim(opt$tsv_name, sep = "	", header = TRUE, stringsAsFactors = FALSE)

if (opt$analyte_name == "all") {
  analyte_names <- df$short_name
} else {
  analyte_names <- c(opt$analyte_name)
}

print(paste("Loading mzML file:", opt$file_name))
ms_data <- readMSData(opt$file_name, mode = "onDisk")
print("File loaded!")

all_results <- list()

for (analyte in analyte_names) {
  print(paste("Processing analyte:", analyte))
  analyte_data <- get_analyte_data(df, analyte)
  mz1 <- assign_mz1(analyte_data, opt$msLevel)
  mz_tol <- (opt$mz_tol_ppm * mz1) / 10^6
  rt_target <- analyte_data$rt_target

  xic <- compute_xic(ms_data, mz1, rt_target, opt$rt_tol_sec, mz_tol, opt$msLevel)
  rt_values <- rtime(xic[[1]])
  int <- intensity(xic[[1]])

  valid_int_idx <- which(!is.na(int) & int > 0)
  if (length(valid_int_idx) > 0) {
    max_idx <- valid_int_idx[which.max(int[valid_int_idx])]
    rt_obs <- rt_values[max_idx]
    intensity_mz1 <- int[max_idx]
  } else {
    rt_obs <- NA
    intensity_mz1 <- NA
  }

  ppp_value <- calculate_ppp_custom(rt_values, int, threshold_ratio = 0.01)
  fwhm_value <- calculate_fwhm(rt_values, int)
  area_mz1 <- compute_xic_area(rt_values, int)
  mz_exp_result <- find_mz_for_intensity_obs(ms_data, c(mz1 - mz_tol, mz1 + mz_tol), rt_obs, intensity_mz1)

  sample_name <- tools::file_path_sans_ext(basename(opt$file_name))
  result_row <- function_create_output_df(
    sample_name, analyte, rt_target, rt_obs, mz1, mz_exp_result,
    intensity_mz1, area_mz1, fwhm_value, ppp_value
  )

  all_results[[analyte]] <- result_row

  cat(paste0("Result for '", analyte, "' added to cumulative TSV output.
"))
}

output_df <- do.call(rbind, all_results)

# Get the retention time and intensity values
rt_values <- rtime(xic[[1]])  
int <- intensity(xic[[1]])  
   
if (opt$msLevel == 2) {
    print(paste("XIC data in MS2:", min(rt_values)))    
    print(paste("Minimum retention time:", min(rt_values)))
    print(paste("Maximum retention time:", max(rt_values)))
    print(paste("Retention time with highest intensity:", rt_values[which.max(int)]))
    print(paste("Maximum intensity in XIC:", max(int, na.rm = TRUE)))
}

if ((plot_xic_ms1 | plot_xic_ms2) && !is.null(plot_output_path)) {
    plot_xic(rt_values, int, plot_output_path)
}

# Verify if there are valid values
valid_int_idx <- which(!is.na(int) & int > 0)

if (length(valid_int_idx) > 0) { # If there are valid values
  max_idx <- valid_int_idx[which.max(int[valid_int_idx])] # Get the index of the maximum intensity
  rt_obs <- rt_values[max_idx] # Get the retention time of the maximum intensity
  intensity_mz1 <- int[max_idx] # Get the intensity of the maximum intensity
} else {
  rt_obs <- NA
  intensity_mz1 <- NA
}

ppp_value <- calculate_ppp_custom(rt_values, int, threshold_ratio = 0.01)
cat("Points Per Peak (PPP):", ppp_value, "\n")
fwhm_value <- calculate_fwhm(rt_values, int)
cat("FWHM (seconds):", fwhm_value, "\n")

# Compute the area under the curve
valid_idx1 <- which(!is.na(rt_values) & !is.na(int)) # Get the valid indexes
area_mz1 <- compute_xic_area(rt_values, int)
log2_area <- ifelse(!is.na(area_mz1) & area_mz1 > 0, log2(area_mz1), NA) # Compute the log2 of the area
delta_mz <- c(mz1 - mz_tol, mz1 + mz_tol) # Compute the delta m/z
mz_exp_result <- find_mz_for_intensity_obs(ms_data, c(mz1 - mz_tol, mz1 + mz_tol), rt_obs, intensity_mz1)
ppm_error <- ((mz_exp_result - mz1) / mz1) * 1e6 # Compute the ppm error

# Print the results
print(paste("Expected retention time (RT target):", rt_target))
print(paste("Observed retention time:", rt_obs))
print(paste("Difference in RT (dRT):", rt_obs - rt_target))
print(paste("Mass m/z1:", mz1))
print(paste("Maximum intensity (M+0):", format_custom(intensity_mz1)))
print(paste("Area M+0:", format_custom(area_mz1)))
print(paste("Total calculated area:", format_custom(area_mz1)))
print(paste("Log2 of total area:", format_custom(log2_area)))
print(paste("Experimental mass found:", mz_exp_result))
print(paste("Difference in m/z (dmz in Da):", mz_exp_result - mz1))
print(paste("Difference in m/z (dmz in ppm):", ppm_error))

# Extract sample name from the input file path
sample_name <- tools::file_path_sans_ext(basename(opt$file_name))

# Ensure output directory exists
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# Force scientific notation globally
options(scipen = 999)

# Define the correct output file path
output_file <- file.path(opt$output_dir, paste0(sample_name, ".tsv"))

if (file.exists(output_file) && !opt$overwrite_tsv) {
  existing_data <- read.delim(output_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  if (!identical(colnames(existing_data), colnames(output_df))) {
    stop("Column mismatch between existing TSV and new output. Consider overwriting.")
  }
  output_df <- rbind(existing_data, output_df)

  # Ensure all numeric columns stay numeric
  numeric_cols <- c("Expected_RT_sec", "Observed_RT_sec", "dRT_sec", "Expected_MZ",
                    "Observed_MZ", "dmz_ppm", "dmz_Da", "Log2_Total_Area")
  output_df[numeric_cols] <- lapply(output_df[numeric_cols], as.numeric)

  # Apply scientific notation ONLY to Max_Intensity and Area
  output_df$Max_Intensity <- format(output_df$Max_Intensity, scientific = TRUE)
  output_df$Area <- format(output_df$Area, scientific = TRUE)
}

write_results(output_df, opt$file_name)

# Write versions.yml for Nextflow traceability
version_info <- list(
  msnbasexic_script = "msnbasexic.R",
  generated_at = as.character(Sys.time()),
  R_version = paste(R.version$major, R.version$minor, sep = "."),
  MSnbase_version = as.character(packageVersion("MSnbase")),
  ggplot2_version = as.character(packageVersion("ggplot2")),
  pracma_version = as.character(packageVersion("pracma"))
)

# Write key-value pairs in YAML-style format
lines <- unlist(lapply(names(version_info), function(k) paste0(k, ": ", version_info[[k]])))
writeLines(lines, "versions.yml")

print("versions.yml generated.")

print(paste("Updated TSV file:", output_file))