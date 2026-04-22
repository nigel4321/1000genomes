# 1000 Genomes Project — Working Notes

## Project Overview
Bioinformatics/genetics analysis project working with **1000 Genomes Project Phase 3** data.

## Local Data Files

All files are in `/home/nigel/1000genomes/`.

| File | Size | Description |
|------|------|-------------|
| `ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | ~329 MB | Chr19 genotypes |
| `ALL.chr19...vcf.gz.tbi` | ~54 KB | Tabix index for chr19 |
| `ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | ~312 MB | Chr20 genotypes |
| `ALL.chr20...vcf.gz.tbi` | ~55 KB | Tabix index for chr20 |
| `ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | ~200 MB | Chr21 genotypes |
| `ALL.chr21...vcf.gz.tbi` | ~35 KB | Tabix index for chr21 |
| `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | ~196 MB | Chr22 genotypes |
| `ALL.chr22...vcf.gz.tbi` | ~35 KB | Tabix index for chr22 |

## File Format Details

- **Format**: VCF (Variant Call Format), bgzipped (`.vcf.gz`)
- **Index**: Tabix (`.tbi`) — enables fast random access by genomic region
- **Pipeline**: `phase3_shapeit2_mvncall_integrated_v5b` — SHAPEIT2 phasing + MVNcall integration, v5b
- **Date stamp**: `20130502` (1000 Genomes Phase 3 freeze date)
- **Samples**: ALL superpopulation (all 2,504 samples across 26 populations)
- **Genome build**: GRCh37 / hg19

## Key Concepts

- **VCF**: stores SNPs, indels, and structural variants with per-sample genotypes
- **bgzip + tabix**: compressed format that allows region queries without decompressing the whole file
- **Phase 3**: final 1000 Genomes release — 2,504 individuals, 84.4M variants

## Common Tools for This Data

```bash
# Inspect header
bcftools view -h ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | head -50

# Query a region
bcftools view -r 22:1000000-2000000 ALL.chr22...vcf.gz

# Count variants
bcftools stats ALL.chr22...vcf.gz | grep "^SN"

# List samples
bcftools query -l ALL.chr22...vcf.gz

# Extract specific fields
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ALL.chr22...vcf.gz

# Filter by allele frequency
bcftools view -q 0.05:minor ALL.chr22...vcf.gz
```

Other relevant tools: `tabix`, `plink`, `vcftools`, `htslib`, Python `cyvcf2` / `pysam`.

## Notes
<!-- Add project-specific notes, goals, and findings here -->
