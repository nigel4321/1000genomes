# Carrier Report — rs12913832

## Summary

- **Carrier records (rows in `carriers.tsv`):** 1274
- **Unique individuals carrying the alt allele:** 637
- **Heterozygotes (dosage = 1):** 772
- **Homozygous alt (dosage = 2):** 502
- **Files containing the variant:** 2
- **Allele-count integrity check:** het + 2·hom = 772 + 2·502 = **1776**

## Per-file breakdown

| File | Heterozygotes | Homozygous alt | Total carriers |
|---|---|---|---|
| `ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | 386 | 251 | 637 |
| `rs12913832.vcf.gz` | 386 | 251 | 637 |

## Output files

- `carriers.tsv` — one row per (file, sample) where the sample carries the alt allele.
  - Columns: `file, sample, variant_id, chrom, pos, ref, alt, genotype, alt_dosage`

## First 10 carriers (preview)

| sample | genotype | alt_dosage |
|---|---|---|
| HG00096 | 1|1 | 2 |
| HG00097 | 1|1 | 2 |
| HG00099 | 1|1 | 2 |
| HG00100 | 1|1 | 2 |
| HG00101 | 1|1 | 2 |
| HG00102 | 1|1 | 2 |
| HG00103 | 1|1 | 2 |
| HG00105 | 1|1 | 2 |
| HG00106 | 1|1 | 2 |
| HG00107 | 0|1 | 1 |

