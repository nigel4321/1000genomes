#!/usr/bin/env bash
# Extract rs12913832 (HERC2, eye-color SNP) from 1000 Genomes Phase 3 chr15.
# Build: GRCh37/hg19 — position 15:28365618, ref A, alt G.

set -euo pipefail

VCF="${1:-ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz}"
RSID="rs12913832"
REGION="15:28365618-28365618"
OUT_PREFIX="${RSID}"

if ! command -v bcftools >/dev/null; then
    echo "bcftools not found. Install with: sudo apt install bcftools" >&2
    exit 1
fi

if [[ ! -f "$VCF" ]]; then
    echo "VCF not found: $VCF" >&2
    echo "Download from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/" >&2
    exit 1
fi

# 1. Fast region query via tabix index, then filter by REF/ALT.
#    The 1000G Phase 3 VCFs were frozen against a 2013 dbSNP snapshot and
#    may not carry the rs12913832 ID, so we identify the variant by its
#    alleles (A>G) rather than by rsID.
REF_ALLELE="A"
ALT_ALLELE="G"
bcftools view -r "$REGION" -i "REF=\"$REF_ALLELE\" & ALT=\"$ALT_ALLELE\"" \
    "$VCF" -Oz -o "${OUT_PREFIX}.vcf.gz"
bcftools index -t "${OUT_PREFIX}.vcf.gz"

echo "--- Variant record ---"
bcftools view -H "${OUT_PREFIX}.vcf.gz" | cut -f1-8

# 2. Per-sample genotype table (sample, GT, dosage of alt allele)
echo "--- Writing genotype table: ${OUT_PREFIX}.genotypes.tsv ---"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GT]\n' \
    "${OUT_PREFIX}.vcf.gz" \
    | tr '\t' '\n' | grep '=' \
    | awk -F= 'BEGIN{OFS="\t"; print "sample","GT","alt_dosage"}
               {
                 gt=$2; n=0
                 if (gt ~ /1/) { split(gt, a, /[|\/]/); for (i in a) if (a[i]=="1") n++ }
                 print $1, gt, n
               }' > "${OUT_PREFIX}.genotypes.tsv"

# 3. Summary: allele frequency across the full 2,504-sample cohort
echo "--- Allele frequency ---"
bcftools +fill-tags "${OUT_PREFIX}.vcf.gz" -- -t AF,AC,AN \
    | bcftools query -f 'AC=%INFO/AC\tAN=%INFO/AN\tAF=%INFO/AF\n'

echo "Done. Outputs:"
echo "  ${OUT_PREFIX}.vcf.gz       — single-variant VCF"
echo "  ${OUT_PREFIX}.vcf.gz.tbi   — index"
echo "  ${OUT_PREFIX}.genotypes.tsv — per-sample genotypes"
