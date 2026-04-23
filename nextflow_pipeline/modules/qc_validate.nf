process QC_VALIDATE {
    tag "$name"

    input:
    tuple val(name), path(vcf), path(tbi)

    output:
    path "${name}.qc.json", emit: json
    tuple val(name), path(vcf), path(tbi), emit: validated

    script:
    def strict = params.strict_qc ? '--strict' : ''
    """
    qc_validate.py \\
        --vcf  "${vcf}" \\
        --name "${name}" \\
        --out  "${name}.qc.json" \\
        ${strict}
    """
}
