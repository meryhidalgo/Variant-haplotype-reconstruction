# Variant-haplotype-reconstruction

Pipeline for reconstructing haplotypes associated with genomic variants using phased long-read sequencing data.

This repository contains Bash and Python scripts to:

- extract phased SNVs from WhatsHap-phased VCF files
- identify phase blocks associated with a genomic locus of interest
- reconstruct local haplotypes across samples
- compare cis-associated variants between individuals
- investigate shared haplotypes and potential founder effects

The workflow is compatible with multiple variant types, including:

- single nucleotide variants (SNVs)
- small insertions/deletions (indels)
- repeat expansions (STRs)
- structural variant-associated loci

The pipeline was designed for phased long-read sequencing datasets generated with technologies such as ONT or PacBio.

---

# Overview

The workflow reconstructs haplotypes surrounding a genomic position of interest by leveraging phased variants contained within the same phase block.

The strategy consists of:

1. Identifying the phase set (PS) associated with a target locus
2. Extracting phased heterozygous variants linked to that phase block
3. Defining shared haplotypic intervals across samples
4. Comparing phased variants among individuals
5. Identifying shared cis-associated alleles and candidate founder haplotypes

This enables comparison of haplotype backgrounds associated with recurrent pathogenic variants or genomic rearrangements.

---

# Requirements

## Software

- Bash
- Python >= 3.8
- bcftools
- htslib
- tabix
- jq

## Python dependencies

- pandas

## Suggested conda environment

```bash
conda create -n haplotyping \
    -c conda-forge \
    -c bioconda \
    bcftools htslib tabix jq pandas
```

---

# Input files

## 1. Phased VCFs directory

Compressed and indexed phased VCF files generated with tools such as WhatsHap. They are expected to be named as: 

```text
sample_whatshap.vcf.gz
sample_whatshap.vcf.gz.tbi
```

## 2. Variant annotation file

Name of 3 


## 3. Genomic position of interest

Genomic coordinate associated with the target variant.

Example:

```bash
POS="chr16:89511445"
```

This position may correspond to SNV, indel, STR expansion etc, but must be reported in the vcf files. 

## 4. Variant annotation file

Annotation file containing rsIDs or genomic annotations.

Example:

```text
hg38_891-092M.summary.tsv
```

---

# Workflow

## Step 1 — Extract phase block-associated heterozygous variants

For each sample:

- identify the phase set (PS) linked to the locus of interest
- extract phased heterozygous variants belonging to the same haplotype block

Output:

```text
haplotyping/txtfiles_heteros/
```

---

## Step 2 — Define shared haplotypic interval

Script:

```bash
python3 5-1-define_block.py
```

This step:

- identifies overlapping phased regions across samples
- defines the shared haplotypic interval
- retrieves cis-associated heterozygous variants

Returned values include:

- interval start
- interval end
- phased heterozygous variant set

---

## Step 3 — Extract phased variants within the interval

All phased variants located within the shared interval are extracted for each sample.

Output:

```text
haplotyping/txtfiles_all/
```

---

## Step 4 — Compare haplotypes across samples

Script:

```bash
python3 5-2-search_variants.py
```

This step identifies:

- shared phased variants
- common cis-associated alleles
- candidate shared haplotypes
- potential founder backgrounds

---

# Example execution

```bash
bash run_haplotyping.sh
```

---

# Applications

- Founder effect analysis
- Haplotype reconstruction
- Long-read phasing studies
- Cis-variant identification
- Shared ancestry analysis
- Variant background characterization

---

# Notes

- The workflow assumes phased genotypes (`0|1`, `1|0`)
- Phase set (`PS`) information is required
- Long-read sequencing data is strongly recommended
- Compatible with any genomic locus represented within phased blocks

---

# Repository structure

```text
variant-haplotype-reconstruction/
│
├── scripts/
│   ├── run_haplotyping.sh
│   ├── 5-1-define_block.py
│   └── 5-2-search_variants.py
│
├── example/
│
├── docs/
│
├── environment.yml
│
└── README.md
```

---

# Future improvements

Potential extensions include:

- haplotype visualization
- BAM haplotag integration
- multi-sample VCF support
- structural variant integration
- graphical haplotype block representation
- automated founder haplotype detection

---

# Citation

If you use this repository, please cite the corresponding study or repository release.

## 3. Variant annotation file

Annotation file containing rsIDs or genomic annotations.

Example:
