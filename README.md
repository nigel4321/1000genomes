# 1000 Genomes — rs12913832 Extraction

## Goal

Extract the **rs12913832** variant from the 1000 Genomes Project Phase 3 data
and produce a per-sample genotype table suitable for downstream analysis
(e.g. population-stratified allele-frequency work).

rs12913832 is a single-nucleotide polymorphism in intron 86 of the **HERC2**
gene, upstream of **OCA2**. It is the primary genetic determinant of
blue-vs-brown eye colour in Europeans — the `G` allele is associated with
blue eyes.

| Field | Value |
|---|---|
| rsID | rs12913832 |
| Gene | HERC2 (intron 86) |
| Build | GRCh37 / hg19 |
| Location | `chr15:28,365,618` |
| Alleles | A (ref) / G (alt) |
| Global AF (G) | ≈ 0.28 |
| EUR AF (G) | ≈ 0.78 |
| AFR AF (G) | ≈ 0.02 |

## Data

The variant lives on chromosome 15, so the `ALL.chr19…` – `ALL.chr22…` files
already in this directory do **not** contain it. Download chr15 from the
1000 Genomes FTP:

```bash
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
FILE="ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

wget "$BASE/$FILE" "$BASE/$FILE.tbi"
```

All autosomes in the 20130502 release use the `v5b` suffix. `.vcf.gz` and
`.tbi` files are covered by `.gitignore` and will not be committed.

## Variant browsers

Useful for sanity-checking coordinates, alleles, and flanking context:

| Resource | URL |
|---|---|
| dbSNP | https://www.ncbi.nlm.nih.gov/snp/rs12913832 |
| Ensembl (GRCh37) | https://grch37.ensembl.org/Homo_sapiens/Variation/Explore?v=rs12913832 |
| UCSC (hg19) | https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=rs12913832 |

Programmatic lookup via NCBI E-utilities:

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=snp&id=12913832&retmode=json" \
  | jq '.result."12913832" | {chr, pos: .chrpos, pos37: .chrpos_prev_assm}'
```

## Script: `extract_rs12913832.sh`

Usage:

```bash
./extract_rs12913832.sh                  # uses default chr15 VCF in cwd
./extract_rs12913832.sh path/to/vcf.gz   # explicit path
```

What it does:

1. **Region query** via tabix (`bcftools view -r 15:28365618-28365618`) — fast
   random access, no full-file scan.
2. **rsID filter** (`-i 'ID="rs12913832"'`) — guards against other variants
   at the same position.
3. Writes a single-variant VCF `rs12913832.vcf.gz` + `.tbi` index.
4. Writes a per-sample genotype table `rs12913832.genotypes.tsv` with
   columns: `sample`, `GT` (e.g. `0|0`, `0|1`, `1|1`), `alt_dosage` (0/1/2).
5. Prints global allele frequency via `bcftools +fill-tags`.

Requires `bcftools` (tested with 1.19).

## Troubleshooting

### Symptom: `rs12913832.genotypes.tsv` has a header but no rows

The header is emitted by awk's `BEGIN{}` regardless of input — an empty TSV
means `bcftools query` returned zero rows. Work through these checks:

```bash
VCF=ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# 1. Is anything at the expected position?
bcftools view -r 15:28365618 -H "$VCF" | head

# 2. Widen the window — 1 kb around the target
bcftools view -r 15:28365000-28366000 -H "$VCF" | awk '{print $1,$2,$3,$4,$5}'

# 3. Authoritative full-file search (slow)
bcftools view "$VCF" | grep -m1 rs12913832

# 4. Confirm the file really is chr15 and check contig naming (15 vs chr15)
bcftools view -h "$VCF" | grep "^##contig" | head
bcftools index -s "$VCF"
```

Interpreting results:

- **Step 4 shows contigs other than `15`** → wrong file. 1000G Phase 3 uses
  unprefixed contig names (`15`, not `chr15`). Re-download chr15.
- **Step 1 is empty but step 2 shows neighbouring variants** → the variant
  may have been repositioned in a later dbSNP build; update `REGION` in the
  script to the coordinate reported in step 2.
- **Step 1 shows the variant but the ID field is e.g. `rs12913832;rsXXXX`**
  → strict `ID="rs12913832"` match fails. Change the filter in the script
  to a regex match:
  ```bash
  bcftools view -r "$REGION" -i 'ID~"rs12913832"' "$VCF"
  ```
- **Step 3 finds the variant but steps 1 and 2 are empty** → the tabix
  index is stale or mismatched. Re-download the `.tbi` file or rebuild with
  `tabix -p vcf "$VCF"`.

## Follow-on analysis

To stratify genotype counts by super-population, join
`rs12913832.genotypes.tsv` against the 1000G sample panel on the `sample`
column:

```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

The panel file maps each sample ID to population (`pop`) and
super-population (`super_pop`: AFR, AMR, EAS, EUR, SAS).
