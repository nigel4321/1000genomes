process METADATA_REPORT {
    publishDir params.outdir, mode: 'copy'

    input:
    path metas

    output:
    path "metadata_report.md"

    script:
    """
    mkdir -p _metas
    for f in ${metas}; do cp "\$f" _metas/; done
    build_metadata_report.py --input-dir _metas --output metadata_report.md
    """
}

process CARRIER_REPORT {
    publishDir params.outdir, mode: 'copy'

    input:
    path carrier_tsvs

    output:
    path "carriers.tsv"
    path "carriers_report.md"

    script:
    """
    mkdir -p _carriers
    for f in ${carrier_tsvs}; do cp "\$f" _carriers/; done
    build_carrier_report.py \\
        --input-dir    _carriers \\
        --output-tsv   carriers.tsv \\
        --output-md    carriers_report.md \\
        --variant-name "${params.variant_name}"
    """
}

process VARIANT_REPORT {
    publishDir params.outdir, mode: 'copy'

    input:
    path results

    output:
    path "variant_report.md"

    script:
    """
    mkdir -p _results
    for f in ${results}; do cp "\$f" _results/; done
    build_variant_report.py \\
        --input-dir     _results \\
        --output        variant_report.md \\
        --variant-name  "${params.variant_name}" \\
        --variant-chrom "${params.variant_chrom}" \\
        --variant-pos   ${params.variant_pos} \\
        --variant-ref   "${params.variant_ref}" \\
        --variant-alt   "${params.variant_alt}" \\
        --min-af        ${params.variant_min_af} \\
        --max-af        ${params.variant_max_af}
    """
}
