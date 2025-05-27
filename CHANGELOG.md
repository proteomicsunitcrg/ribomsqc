# nf-core/ribomsqc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - 2025-03-26

### `Added`

* [a74c9b1](https://github.com/proteomicsunitcrg/ribomsqc/commit/a74c9b1) - Added new quantification markers to XIC line plots for enhanced interpretability
* [b2d9f6e](https://github.com/proteomicsunitcrg/ribomsqc/commit/b2d9f6e) - Accumulate TSV export with overwrite toggle
* [bd78c6a](https://github.com/proteomicsunitcrg/ribomsqc/commit/bd78c6a) - Add MS level to section\_name in MultiQC JSON output
* [4a1f2f3](https://github.com/proteomicsunitcrg/ribomsqc/commit/4a1f2f3) - Enhance XIC plots with FWHM and PPP threshold visualization
* [e4b8f32](https://github.com/proteomicsunitcrg/ribomsqc/commit/e4b8f32) - Add support for --analyte\_name all to process all TSV entries
* [5563206](https://github.com/proteomicsunitcrg/ribomsqc/commit/55632066ab4eb1191a52b71e887d979e1e8ca6e2) - Support for `MultiQC` as module to aggregate results.
* [9a16ef6](https://github.com/proteomicsunitcrg/ribomsqc/commit/9a16ef67a1a4b3d675ddbfc5615158b321eba7bc) - Custom `MSNBASEXIC` module to extract and plot XICs.
* [b7ea318](https://github.com/proteomicsunitcrg/ribomsqc/commit/b7ea318d7ac3422e5bf6e87a769af6284d3024aa) - Support for ThermoRawFileParser module to convert RAW to mzML.
* [8d4c4e3](https://github.com/proteomicsunitcrg/ribomsqc/commit/8d4c4e368e0f774b246002e5e73c8bf53aab4391) - Initial implementation of the pipeline.

### `Changed`

* [c35b1fe](https://github.com/proteomicsunitcrg/ribomsqc/commit/c35b1fe) - Improved logging clarity and enhanced plot aesthetics
* [5c9eabc](https://github.com/proteomicsunitcrg/ribomsqc/commit/5c9eabc) - Corrected MS2 peak extraction logic with accurate precursor matching
* [517fa6e](https://github.com/proteomicsunitcrg/ribomsqc/commit/517fa6e) - Refined internal logic of the XIC extraction algorithm

### `Fixed`

* [20fa8c7](https://github.com/proteomicsunitcrg/ribomsqc/commit/20fa8c7) - Resolved issue causing missing JSON output entries for certain analytes
* [d8f4c5b](https://github.com/proteomicsunitcrg/ribomsqc/commit/d8f4c5b) - Harmonized XIC plot styling, layout and retention time interval handling
* [e3ac9e2](https://github.com/proteomicsunitcrg/ribomsqc/commit/e3ac9e2) - Preserve MS level in mergejsons
* [f6c8e1d](https://github.com/proteomicsunitcrg/ribomsqc/commit/f6c8e1d) - Fixed FWHM interpolation logic

### `Dependencies`

| Dependency                 | Version |
| -------------------------- | ------- |
| `ThermoRawFileParser`      | 1.4.5   |
| `bioconductor-msnbase (R)` | 2.32.0  |
| `MultiQC`                  | 1.28    |
