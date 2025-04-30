#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(jsonlite)
  library(tools)
})

# Define input folder for staged JSONs
input_dir <- "merge_jsons_input"

# Validate that the input directory exists
if (!dir.exists(input_dir)) {
  stop("Input directory 'merge_jsons_input' does not exist.")
}

# List all candidate *_mqc.json files inside the input directory
json_files <- list.files(path = input_dir, pattern = "*_mqc.json$", full.names = TRUE)

# Abort if no JSONs are found
if (length(json_files) == 0) {
  stop("No JSON files found in input directory.")
}

cat("[INFO] Found JSON input files:\n")
print(json_files)

# Try to read the metric IDs from each file (with error handling)
extract_metric_id <- function(f) {
  tryCatch({
    fromJSON(f)$id
  }, error = function(e) {
    warning(paste("Skipping invalid JSON file:", f, "->", e$message))
    return(NA)
  })
}

metric_ids <- unique(na.omit(sapply(json_files, extract_metric_id)))

# Abort if no valid metric IDs were parsed
if (length(metric_ids) == 0) {
  stop("No valid JSON files with a 'id' field found.")
}

# Check if all expected merged outputs already exist
existing_outputs <- list.files(pattern = "_merged_mqc\\.json$")

expected_outputs <- paste0(metric_ids, "_merged_mqc.json")

if (all(expected_outputs %in% existing_outputs)) {
  cat("[INFO] All merged JSONs already exist. Skipping merge.\n")
  quit(status = 0)
}

# If some outputs exist and others do not, raise an error
if (any(expected_outputs %in% existing_outputs)) {
  stop("Some merged JSONs already exist. Please remove or clean the output directory before reprocessing.")
}

# Initialize the merged data structure
merged_metrics <- list()

# Process each JSON file and group values by metric → analyte → sample
for (file in json_files) {
  cat("[INFO] Processing file:", file, "\n")
  data <- tryCatch({
    fromJSON(file)
  }, error = function(e) {
    warning(paste("Skipping invalid JSON:", file, "->", e$message))
    return(NULL)
  })

  if (is.null(data) || is.null(data$id) || is.null(data$data)) {
    warning(paste("Skipping malformed JSON (missing 'id' or 'data'):", file))
    next
  }

  metric <- data$id

  if (!metric %in% names(merged_metrics)) {
    merged_metrics[[metric]] <- list()
  }

  for (analyte in names(data$data)) {
    if (!analyte %in% names(merged_metrics[[metric]])) {
      merged_metrics[[metric]][[analyte]] <- list()
    }

    for (sample in names(data$data[[analyte]])) {
      merged_metrics[[metric]][[analyte]][[sample]] <- data$data[[analyte]][[sample]]
    }
  }
}

# Write merged JSONs, one per metric
for (metric in names(merged_metrics)) {
  merged_json <- list(
    id = metric,
    section_name = metric,
    description = paste(metric, "values across samples"),
    plot_type = "linegraph",
    pconfig = list(
      id = paste0(metric, "_plot"),
      title = metric,
      xlab = "Sample",
      ylab = paste0(metric, " (unit)"),
      xlab_format = "category",
      showlegend = TRUE
    ),
    data = merged_metrics[[metric]]
  )

  out_file <- paste0(metric, "_merged_mqc.json")
  cat("[INFO] Writing merged JSON:", out_file, "\n")
  write_json(merged_json, path = out_file, auto_unbox = TRUE, pretty = TRUE)
}