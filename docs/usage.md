# nf-core/ribomsqc: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/ribomsqc/usage](https://nf-co.re/ribomsqc/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

`nf-core/ribomsqc` is a quality control (QC) pipeline designed for monitoring mass spectrometry performance in ribonucleoside analysis. It parses RAW files, performs XIC extraction using analyte definitions, and generates summary plots with MultiQC.

## Samplesheet input

Before running the pipeline, you need to provide a samplesheet CSV file containing metadata about the files to process.

```bash
--input '[path to samplesheet.csv]'
```

The samplesheet must contain a header and two columns:

```csv title="samplesheet.csv"
id,raw_file
20240315_QCN1_001_03_New_STD,/full/path/to/file.raw
```

- `id`: Sample identifier
- `file`: Path to the RAW file

## Running the pipeline

Typical command:

```bash
nextflow run nf-core/ribomsqc \
  --input /home/proteomics/mydata/csv/samplesheet.csv \
  --analytes_tsv /home/proteomics/mydata/tsv/qcn1.tsv \
  --analyte m3C \
  --rt_tolerance 150 \
  --mz_tolerance 20 \
  --ms_level 2 \
  --plot_xic_ms1 false \
  --plot_xic_ms2 false \
  --plot_output_path xic_plot \
  --overwrite_tsv true \
  --outdir results \
  -profile singularity
```

Alternatively to specifying parameters on the command line, you can define them in a separate parameters file (e.g., `params.json` or `params.yaml`). This can help simplify execution and improve reproducibility, especially for complex or repeated runs.

## Minimal params file example

```yaml title="params.yaml"
input: /home/proteomics/mydata/csv/samplesheet.csv
analytes_tsv: /home/proteomics/mydata/tsv/qcn1.tsv
analyte: m3C
rt_tolerance: 150
mz_tolerance: 20
ms_level: 2
plot_xic_ms1: false
plot_xic_ms2: false
plot_output_path: xic_plot
overwrite_tsv: true
outdir: results
```

Then run with:

```bash
nextflow run nf-core/ribomsqc -profile singularity -params-file params.yaml
```

## MultiQC Integration

The pipeline integrates a custom subworkflow to run [MultiQC](https://multiqc.info/). It collects `.tsv` files generated by the `MSNBASEXIC` module and summarises them. Output is stored in `${params.outdir}/multiqc/`.

## Reproducibility

Use `-r` to specify a pipeline version:

```bash
nextflow run nf-core/ribomsqc -r 1.0.1dev ...
```

Update pipeline:

```bash
nextflow pull nf-core/ribomsqc
```

## Core Nextflow options

- `-profile docker|singularity|conda|podman|...`: Choose environment
- `-resume`: Resume from previous run
- `-params-file`: Load parameters from YAML or JSON file
- `-c`: Load additional config for cluster resources, etc.

## Tips

- Use `-profile singularity` for reproducibility
- The pipeline is compatible with [Wave containers](https://seqera.io/wave/) for dynamic container resolution
- For cluster configs, use `nf-core/configs`

## Example output

```bash
results/
├── multiqc/
│   └── multiqc_report.html
├── xic_results/
│   ├── sample1.tsv
│   └── ...
└── pipeline_info/
    └── nf_core_ribomsqc_software_versions.yml
```

---

Generated with :heart: by [nf-core](https://nf-co.re) and adapted for custom analyte QC workflows.