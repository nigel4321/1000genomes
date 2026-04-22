process INSPECT_VCF {
    tag "$name"

    input:
    tuple val(name), path(vcf), path(tbi)

    output:
    path "${name}.meta.json"

    script:
    """
    inspect_vcf.py "${vcf}" "${name}" > "${name}.meta.json"
    """
}
