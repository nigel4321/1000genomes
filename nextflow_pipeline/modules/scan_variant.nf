process SCAN_VARIANT {
    tag "$name"

    input:
    tuple val(name), path(vcf), path(tbi)

    output:
    path "${name}.variant.json", emit: json
    path "${name}.carriers.tsv", emit: carriers

    script:
    """
    scan_variant.py \\
        --vcf        "${vcf}" \\
        --name       "${name}" \\
        --out-prefix "${name}" \\
        --chrom      "${params.variant_chrom}" \\
        --pos        ${params.variant_pos} \\
        --ref        "${params.variant_ref}" \\
        --alt        "${params.variant_alt}" \\
        --min-af     ${params.variant_min_af} \\
        --max-af     ${params.variant_max_af}
    """
}
