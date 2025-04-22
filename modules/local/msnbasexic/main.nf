process MSNBASEXIC {
    label 'process_single'
    tag "$analyte_id"

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/bioconductor-msnbase_r-ggplot2_r-optparse_r-pracma_r-readr:83cd263d3bfd0c9e'

input:
    val analyte_id
    val row
    path mzml_file
    path tsv_file

    output:
        tuple val(analyte_id), path("*.tsv"), emit: xic_output
        path "versions.yml", emit: versions

    script:
    def outdir   = "${params.outdir}/xic_results"
    def plot_dir = params.plot_output_path

    """
    msnbasexic.R \\
    --file_name ${mzml_file} \\
    --tsv_name ${tsv_file} \\
    --analyte_name ${analyte_id} \\
    --output_dir ${outdir} \\
    --rt_tol_sec ${params.rt_tolerance} \\
    --mz_tol_ppm ${params.mz_tolerance} \\
    --msLevel ${params.ms_level} \\
    --plot_xic_ms1 ${params.plot_xic_ms1} \\
    --plot_xic_ms2 ${params.plot_xic_ms2} \\
    --plot_output_path ${plot_dir} \\
    --overwrite_tsv ${params.overwrite_tsv}
    """



}