# nf-core/ribomsqc: Output

## Introduction

This document describes the output produced by the pipeline.
After execution, the following directories will be created under the top-level `results/` directory, each containing key output files necessary for downstream QC and reporting.

## Pipeline overview

The pipeline is implemented with [Nextflow](https://www.nextflow.io/) and executes the following major steps:

1. **ThermoRawFileParser** – Conversion of raw instrument files
2. **MSNBASEXIC** – Extract chromatographic peaks and calculate per-analyte QC metrics
3. **Merge JSONs** – Consolidation of individual QC metric JSONs into merged summaries
4. **MultiQC** – Aggregation and visualization of QC metrics
5. **Pipeline information** – Generation of execution reports and parameter logs

---

## ThermoRawFileParser

<details markdown="1">
<summary>Output files</summary>

* `thermorawfileparser/`

  * Converted `.raw` files to `.mzML` format suitable for downstream processing.

</details>

---

## XIC results

<details markdown="1">
<summary>Output files</summary>

* `msnbasexic/`

  * One `.json` file is generated per sample, containing quantitative chromatographic data and quality control (QC) metrics derived from the XIC-based algorithm.

</details>

---

## Merge JSONs

<details markdown="1">
<summary>Output files</summary>

* `mergejsons/`

  * `dmz_ppm_merged_mqc.json`
  * `FWHM_merged_mqc.json`
  * `Log2_Total_Area_merged_mqc.json`
  * `Observed_RT_sec_merged_mqc.json`

Each of these **merged** JSON files consolidates multiple per-sample (or per-file) JSON outputs for the same QC metric into a single array. For example, all individual `dmz_ppm.json` fragments from different samples are combined into `dmz_ppm_merged_mqc.json`. This ensures that each parameter’s values across all samples are gathered in one file, formatted for seamless ingestion by MultiQC and unified visualization of QC trends.

</details>

---

## MultiQC report

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`

  * `multiqc_report.html`: Interactive HTML report aggregating XIC-derived tables and merged JSON metrics.
  * `multiqc_data/`: Directory containing raw data inputs for MultiQC, including the merged QC JSON files.

</details>

---

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`

  * Nextflow execution artifacts: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt`, `pipeline_dag.dot`/`pipeline_dag.svg`
  * Pipeline software versions: `nf_core_ribomsqc_software_versions.yml`
  * Run parameters: `params.json` (only if `-saveParams true` is set)

</details>

For more details on Nextflow execution reports and troubleshooting, refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/tracing.html).
