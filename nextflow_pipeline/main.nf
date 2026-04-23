#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QC_VALIDATE }                     from './modules/qc_validate.nf'
include { INSPECT_VCF }                     from './modules/inspect_vcf.nf'
include { SCAN_VARIANT }                    from './modules/scan_variant.nf'
include { QC_REPORT; METADATA_REPORT; VARIANT_REPORT; CARRIER_REPORT } from './modules/reports.nf'

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

    QC_VALIDATE(vcf_ch)

    // Gate downstream stages on QC completion. In strict mode qc_validate.py
    // exits non-zero on hard failure, so this tuple is never emitted for a
    // broken file. In non-strict mode everything flows through regardless.
    validated_ch = QC_VALIDATE.out.validated

    INSPECT_VCF(validated_ch)
    SCAN_VARIANT(validated_ch)

    QC_REPORT(QC_VALIDATE.out.json.collect())
    METADATA_REPORT(INSPECT_VCF.out.collect())
    VARIANT_REPORT(SCAN_VARIANT.out.json.collect())
    CARRIER_REPORT(SCAN_VARIANT.out.carriers.collect())
}

workflow.onComplete {
    log.info "Pipeline complete — reports in ${params.outdir}"
}
