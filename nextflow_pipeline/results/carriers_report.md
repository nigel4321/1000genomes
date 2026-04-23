# Carrier Report — rs12913832

## Summary

- **Carrier records (rows in `carriers.tsv`):** 4
- **Unique individuals carrying the alt allele:** 4
- **Heterozygotes (dosage = 1):** 1
- **Homozygous alt (dosage = 2):** 3
- **Files containing the variant:** 1
- **Allele-count integrity check:** het + 2·hom = 1 + 2·3 = **7**

## Per-file breakdown

| File | Heterozygotes | Homozygous alt | Total carriers |
|---|---|---|---|
| `ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | 1 | 3 | 4 |

## Output files

- `carriers.tsv` — one row per (file, sample) where the sample carries the alt allele.
  - Columns: `file, sample, variant_id, chrom, pos, ref, alt, genotype, alt_dosage`

## First 10 carriers (preview)

| sample | genotype | alt_dosage |
|---|---|---|
| HG00096 | 1|1 | 2 |
| HG00097 | 1|1 | 2 |
| HG00100 | 1|1 | 2 |
| HG01112 | 0|1 | 1 |

