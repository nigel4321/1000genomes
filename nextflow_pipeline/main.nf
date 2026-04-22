#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INSPECT_VCF }                     from './modules/inspect_vcf.nf'
include { SCAN_VARIANT }                    from './modules/scan_variant.nf'
include { METADATA_REPORT; VARIANT_REPORT; CARRIER_REPORT } from './modules/reports.nf'

def requiredParams = [
    'input',
    'variant_name', 'variant_chrom', 'variant_pos',
    'variant_ref',  'variant_alt',
]
requiredParams.each { k ->
    if (params[k] == null) {
        error "Missing required parameter: --${k} (see README.md)"
    }
}

workflow {
    vcf_ch = Channel.fromPath(params.input, checkIfExists: true)
        .map { vcf ->
            def tbi = file("${vcf}.tbi")
            if (!tbi.exists()) {
                error "Missing tabix index for ${vcf} — expected at ${tbi}"
            }
            def name = vcf.name.replaceAll(/\.vcf\.gz$|\.vcf\.bgz$|\.vcf$/, "")
            tuple(name, vcf, tbi)
        }

    INSPECT_VCF(vcf_ch)
    SCAN_VARIANT(vcf_ch)

    METADATA_REPORT(INSPECT_VCF.out.collect())
    VARIANT_REPORT(SCAN_VARIANT.out.json.collect())
    CARRIER_REPORT(SCAN_VARIANT.out.carriers.collect())
}

workflow.onComplete {
    log.info "Pipeline complete — reports in ${params.outdir}"
}
