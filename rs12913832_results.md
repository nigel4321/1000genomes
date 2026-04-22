# rs12913832 Extraction — Results Summary

## Objective

Extract the **rs12913832** variant (HERC2 intron 86, GRCh37 `chr15:28,365,618`,
A>G) from 1000 Genomes Project Phase 3 genotypes and produce a per-sample
genotype table for downstream population-stratified analysis. This SNP is
the primary determinant of blue-vs-brown eye colour in Europeans.

## Data and tools

| Item | Value |
|---|---|
| Source | `ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` |
| Build | GRCh37 / hg19 |
| Samples | 2,504 (all 26 populations) |
| Tools | `bcftools 1.19`, `tabix` |
| Script | `extract_rs12913832.sh` |

The chr15 VCF was downloaded from the 1000 Genomes EBI FTP; the chr19–chr22
files already present in the working directory do not contain this variant.

## Pipeline

1. Region query via tabix index — `bcftools view -r 15:28365618-28365618`
2. Filter to the A>G biallelic SNP — `-i 'REF="A" & ALT="G"'`
3. Write single-variant VCF + tabix index
4. Emit per-sample TSV (`sample`, `GT`, `alt_dosage`) via `bcftools query` + awk
5. Recompute global allele counts with `bcftools +fill-tags`

## Troubleshooting path

The initial filter `-i 'ID="rs12913832"'` produced an empty output. Diagnostic
steps established:

1. A variant **is** present at `15:28365618` — but its ID field is `.`.
2. A full-file search for the string `rs12913832` returned no matches —
   the 1000G Phase 3 release was frozen against a 2013 dbSNP snapshot and
   does not carry this rsID.
3. The REF/ALT at that position (`A`/`G`) matches the canonical alleles for
   rs12913832 from dbSNP/Ensembl.

**Resolution:** switched the filter from rsID-based to position + allele
match. This is the more robust approach generally — it survives dbSNP
re-numbering and is not dependent on which snapshot the VCF was frozen against.

A secondary observation: the same tabix region also returns two copy-number
variants (`<CN2>` and `<CN0>,<CN2>`) at nearby positions (28343591, 28353161).
The REF/ALT filter correctly excludes these.

## Results

### Record extracted

```
15  28365618  .  A  G  100  PASS
AC=888; AF=0.177316; AN=5008; NS=2504; DP=19161
EAS_AF=0.002; AMR_AF=0.2017; AFR_AF=0.028; EUR_AF=0.6362; SAS_AF=0.0706
AA=A; VT=SNP
```

### Integrity checks

| Check | Expected | Actual | Pass |
|---|---|---|---|
| `AN` | 2 × 2504 = 5008 | 5008 | ✅ |
| `AC` vs recomputed | Match | 888 both | ✅ |
| `AF` vs recomputed | Match | 0.177316 both | ✅ |
| Variant type | SNP | `VT=SNP` | ✅ |
| Missingness | None | `AN=5008` full | ✅ |

### Allele frequencies

| Scope | AF (G allele) |
|---|---|
| Global | 0.177 |
| EUR (European) | **0.636** |
| AMR (Admixed American) | 0.202 |
| SAS (South Asian) | 0.071 |
| AFR (African) | 0.028 |
| EAS (East Asian) | 0.002 |

## Interpretation

- **Ancestral allele is A, derived allele is G** (`AA=A`). The blue-eye–associated
  allele is the derived mutation.
- **European enrichment with near-absence elsewhere** is the diagnostic
  ancestry signature for rs12913832, consistent with the derived G allele
  arising and rising to high frequency along the lineage leading to modern
  Europeans, with subsequent partial introgression into SAS and AMR via
  admixture.
- The moderate AMR frequency (0.20) is expected given European ancestry
  contribution to admixed American populations.
- The ~300-fold difference between EUR (0.636) and EAS (0.002) is one of
  the largest population-frequency differentials in the human genome and is
  routinely cited as a textbook example of local positive selection.

## Determination

The variant at `15:28365618` with REF=A / ALT=G in the 1000G Phase 3 chr15
VCF is confirmed to be **rs12913832**, on three independent lines of evidence:

1. **Coordinate match** — position is identical to the dbSNP/Ensembl GRCh37
   record.
2. **Allele match** — REF/ALT are the canonical A/G.
3. **Frequency fingerprint** — the per-population AF pattern (very high EUR,
   near-zero EAS/AFR) matches published figures for rs12913832 and is
   effectively unique at this locus.

The missing rsID in the ID column is a known limitation of the Phase 3
release snapshot and does not affect identity.

## Outputs

| File | Contents |
|---|---|
| `rs12913832.vcf.gz` | Single-variant VCF |
| `rs12913832.vcf.gz.tbi` | Tabix index |
| `rs12913832.genotypes.tsv` | 2,504 rows: `sample`, `GT`, `alt_dosage` |

## Suggested next steps

- Join `rs12913832.genotypes.tsv` with `integrated_call_samples_v3.20130502.ALL.panel`
  on `sample` to reproduce per-population genotype counts from raw data.
- Test Hardy–Weinberg equilibrium within each super-population as a QC step.
- Extend the same pipeline (position + allele filter) to other pigmentation
  variants (e.g. rs1426654 in SLC24A5, rs16891982 in SLC45A2) for a
  multi-locus ancestry/pigmentation panel.
