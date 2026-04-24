This Markdown document serves as a technical specification for a developer or bioinformatics engineer to build a high-fidelity synthetic genetic data generator. It captures the requirements for VCF structure, biological realism, and the specific nuances of a mixed UK-based cohort.

***

# Technical Specification: High-Fidelity Synthetic Genetic Data Generator

## 1. Project Objective
To develop a scalable pipeline for generating synthetic (mock) human genetic data in VCF format. The data must be indistinguishable from real-world sequencing results for the purpose of benchmarking bioinformatics pipelines, testing clinical variant callers, and performing population genetic research on a mixed UK cohort.

---

## 2. File Format Requirements (VCF 4.2/4.3)
All output must strictly adhere to the **GA4GH** standards. A "clean" VCF is insufficient; it must mimic the technical metadata of a real sequencer.

### 2.1 Header Requirements
* **Contig Lines:** Explicitly define every chromosome (1–22, X, Y, MT) with lengths and assembly source (`##contig=<ID=chr1,length=248956422,assembly=GRCh38>`).
* **Metadata Tags:** Define all `INFO` and `FORMAT` tags.
* **Reference:** Data must be aligned to **GRCh38** (Genome Reference Consortium Human Build 38).

### 2.2 Variant Records (The Body)
* **Genotype Fields (`FORMAT`):**
    * `GT`: Phased genotypes (e.g., `0|1`) to support linkage studies.
    * `DP`: Simulated depth of coverage (stochastic distribution, e.g., Poisson).
    * `GQ`: Genotype Quality scores mimicking confidence variance.
    * `AD`: Allelic Depth showing realistic sampling bias (e.g., 45/55 for a heterozygote).
* **Variant Diversity:**
    * Include SNPs, short Indels (left-aligned), and Structural Variants (SVs) using `SVTYPE` tags.
    * Incorporate multi-allelic sites.

---

## 3. Biological Realism & Linkage
The generator must not treat variants as independent events. It must simulate the evolutionary constraints of the human genome.

### 3.1 Linkage Disequilibrium (LD)
* **Requirement:** Variants must be clustered into LD blocks following a recombination map (e.g., deCODE).
* **Implementation:** Use a **Coalescent Simulation** approach (e.g., `msprime`) or a **Hidden Markov Model** (e.g., Li and Stephens model) to ensure $r^2$ and $D'$ values decay realistically with distance.

### 3.2 Mutation Metrics
* **Ti/Tv Ratio:** Maintain a Transition/Transversion ratio of ~2.1 for WGS.
* **Site Frequency Spectrum (SFS):** Most variants must be rare (Singletons), following the power-law distribution seen in gnomAD.
* **Database Grounding:** Inject known variants from **dbSNP**, **ClinVar** (pathogenic markers), and **COSMIC** (if simulating somatic data).

---

## 4. Population Simulation: Mixed UK Cohort
The data must reflect the demographic diversity of a major UK urban center (e.g., London), requiring complex admixture modeling.

### 4.1 Ancestral Components
Simulate a mosaic genome for each individual based on:
1.  **Northern/Western European (EUR)**
2.  **South Asian (SAS)** - Indian, Pakistani, Bangladeshi lineages.
3.  **African/Caribbean (AFR)** - Higher diversity, shorter LD blocks.

### 4.2 Admixture Modeling
* **Local Ancestry:** Implement an "Admixture Pulse" model where recombination events switch between ancestral templates within a single chromosome.
* **Metadata:** Provide a secondary "Truth Map" (BED format) indicating the ancestral origin of every segment for each synthetic individual.

---

## 5. Recommended Algorithmic Stack

| Component | Recommended Tool/Library |
| :--- | :--- |
| **Ancestry & LD** | `msprime` (Python) - Scalable tree sequence recording. |
| **Variant Injection** | `bcftools` or custom script to overlay ClinVar/dbSNP alleles. |
| **Demographics** | `stdpopsim` - Pre-configured human demographic models. |
| **Error Modeling** | `ART` or `SimNGS` - To simulate sequencing machine noise. |

---

## 6. Validation Suite (Acceptance Criteria)
The developer must provide validation plots for every generated batch:

1.  **PCA Plot:** Synthetic samples must cluster or bridge clusters (admixture) correctly against 1000 Genomes reference samples.
2.  **LD Decay Curve:** A plot of $r^2$ vs. physical distance (kb) must match biological expectations for the target population.
3.  **Variant Statistics:** A summary report showing $Ti/Tv$ ratios, Het/Hom ratios, and Allele Frequency distributions.

---

## 7. Delivery Format
1.  **VCF.gz:** Compressed and indexed (`.tbi`).
2.  **Truth Set:** A BED file detailing the "Golden Truth" (actual variants vs. simulated noise).
3.  **Local Ancestry Map:** For admixed individuals.
