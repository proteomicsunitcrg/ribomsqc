process MSNBASEXIC {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    // Singularity (ociAutoPull via Wave)
    'community.wave.seqera.io/library/bioconductor-msnbase_r-ggplot2_r-optparse_r-pracma_r-readr:83cd263d3bfd0c9e' :
    // Docker/Wave
    'community.wave.seqera.io/library/bioconductor-msnbase_r-ggplot2_r-optparse_r-pracma_r-readr:83cd263d3bfd0c9e' }"

    input:
        tuple val(meta), path(mzml_file), path(tsv_file)
        val analyte
        val rt_tol_sec
        val mz_tol_ppm
        val msLevel
        val plot_xic_ms1
        val plot_xic_ms2
        val plot_output_path
        val overwrite_tsv

    output:
        tuple val(meta), path("*.tsv"), emit: xic_output
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${workflow.projectDir}/modules/local/msnbasexic/resources/users/bin/msnbasexic.R \\
    --file_name ${mzml_file} \\
    --tsv_name ${tsv_file} \\
    --output_dir ${params.outdir}/xic_results \\
    --analyte_name ${analyte} \\
    --rt_tol_sec ${rt_tol_sec} \\
    --mz_tol_ppm ${mz_tol_ppm} \\
    --msLevel ${msLevel} \\
    --plot_xic_ms1 ${plot_xic_ms1} \\
    --plot_xic_ms2 ${plot_xic_ms2} \\
    --plot_output_path ${plot_output_path} \\
    --overwrite_tsv ${overwrite_tsv}
    """
}