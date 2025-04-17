<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-ribomsqc_logo_dark.png">
    <img alt="nf-core/ribomsqc" src="docs/images/nf-core-ribomsqc_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/ribomsqc/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/ribomsqc/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/ribomsqc/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/ribomsqc/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/ribomsqc/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/ribomsqc)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23ribomsqc-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/ribomsqc)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/ribomsqc** is a bioinformatics pipeline that processes `RAW` files from mass spectrometry instruments, converts them to `mzML` format, extracts `XIC`s for selected analytes, generates plots, and summarizes outputs via MultiQC. It is tailored for quality control of ribonucleoside analysis.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

First, prepare a samplesheet with your input data like:

```csv title="samplesheet.csv"
id,raw_file
20240315_QCN1_001_03_New_STD,/full/path/to/file.raw
```

- `id`: Sample identifier
- `file`: Path to the `RAW` file

Now, you can run the pipeline using:

```bash
nextflow run nf-core/ribomsqc \
  --input samplesheet.csv \
  --analytes_tsv qcn1.tsv \
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

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option.

For more information, see the [usage docs](https://nf-co.re/ribomsqc/usage) and [parameters](https://nf-co.re/ribomsqc/parameters).

## Pipeline output

See [results page](https://nf-co.re/ribomsqc/results) for example output and [output docs](https://nf-co.re/ribomsqc/output).

## Credits

nf-core/ribomsqc was originally written by Roger Olivella.

## Contributions and Support

For help, visit [Slack #ribomsqc](https://nfcore.slack.com/channels/ribomsqc) or see [contributing guide](.github/CONTRIBUTING.md).

## Citations

See [`CITATIONS.md`](CITATIONS.md) for tool references.

> Ewels PA *et al.* (2020) _The nf-core framework_. Nat Biotechnol. [doi:10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)