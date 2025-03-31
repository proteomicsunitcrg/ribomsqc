# nf-core/ribomsqc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - 2025-03-26

Initial release of nf-core/ribomsqc, created with the [nf-core](https://nf-co.re/) template.

### `Added`
- [8d4c4e3](https://github.com/proteomicsunitcrg/ribomsqc/commit/8d4c4e368e0f774b246002e5e73c8bf53aab4391) - Initial implementation of the pipeline.
- [b7ea318](https://github.com/proteomicsunitcrg/ribomsqc/commit/b7ea318d7ac3422e5bf6e87a769af6284d3024aa) - Support for ThermoRawFileParser module to convert RAW to mzML.
- [9a16ef6](https://github.com/proteomicsunitcrg/ribomsqc/commit/9a16ef67a1a4b3d675ddbfc5615158b321eba7bc) - Custom `MSNBASEXIC` module to extract and plot XICs.
- [5563206](https://github.com/proteomicsunitcrg/ribomsqc/commit/55632066ab4eb1191a52b71e887d979e1e8ca6e2) - Support for `MultiQC` as module to aggregate results.

### `Dependencies`

| Dependency                 | Version     | 
| -----------------------    | ----------- |
| `ThermoRawFileParser`      | 1.4.5       |
| `bioconductor-msnbase (R)` | 2.32.0      | 
| `MultiQC`                  | 1.28     | 
