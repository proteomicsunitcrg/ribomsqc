/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ribomsqc_pipeline'
include { THERMORAWFILEPARSER } from '../modules/nf-core/thermorawfileparser/main'
include { MSNBASEXIC } from '../modules/local/msnbasexic/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RIBOMSQC {

    take:
    input_ch

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run THERMORAWFILEPARSER
    //
    THERMORAWFILEPARSER(
         input_ch
    )
    
    ch_versions = ch_versions.mix(THERMORAWFILEPARSER.out.versions)

    mzml_ch = THERMORAWFILEPARSER.out.spectra.map { meta, mzml_file ->
        tuple(
            meta,                      
            mzml_file,
            file(params.analytes_tsv)
        )
    }

    tsv_ch              = Channel.value(file(params.analytes_tsv))
    analyte_ch          = Channel.value(params.analyte)
    rt_tol_ch           = Channel.value(params.rt_tolerance)
    mz_tol_ch           = Channel.value(params.mz_tolerance)
    ms_level_ch         = Channel.value(params.ms_level)
    plot_xic_ms1_ch     = Channel.value(params.plot_xic_ms1)
    plot_xic_ms2_ch     = Channel.value(params.plot_xic_ms2)
    plot_output_path_ch = Channel.value(params.plot_output_path)
    overwrite_tsv_ch    = Channel.value(params.overwrite_tsv)

    MSNBASEXIC(
        mzml_ch,
        Channel.value(params.analyte),
        Channel.value(params.rt_tolerance),
        Channel.value(params.mz_tolerance),
        Channel.value(params.ms_level),
        Channel.value(params.plot_xic_ms1),
        Channel.value(params.plot_xic_ms2),
        Channel.value(params.plot_output_path),
        Channel.value(params.overwrite_tsv)
    )

    ch_versions = ch_versions.mix(MSNBASEXIC.out.versions)


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'ribomsqc_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions = ch_versions                     
    spectra  = THERMORAWFILEPARSER.out.spectra
    xic_output = MSNBASEXIC.out.xic_output

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/