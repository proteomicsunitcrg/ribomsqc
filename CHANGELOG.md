# nf-core/ribomsqc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.1dev - 2025-03-26

Initial release of nf-core/ribomsqc, created with the [nf-core](https://nf-co.re/) template.

### `Added`
- Initial implementation of the pipeline.
- Support for ThermoRawFileParser module to convert RAW to mzML.
- Custom `MSNBASEXIC` module to extract and plot XICs.
- Parameters for XIC configuration: analyte, mz/rt tolerance, msLevel, plotting.
- Support for `MultiQC` as subworkflow to aggregate results.
- Singularity and Docker support.
- Input schema validation via `nextflow_schema.json`.
- Custom samplesheet CSV input with three required columns.
- Version tracking via `software_versions.yml`.

### `Fixed`
- N/A – First release.

### `Dependencies`
- ThermoRawFileParser
- bioconductor-msnbase (R)
- MultiQC
- nf-core/tools v3.2.0+

### `Deprecated`
- Nothing deprecated in this version.