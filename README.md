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

## 1. Phased VCFs

Compressed and indexed phased VCF files generated with tools such as WhatsHap.

Example:

```text
sample_whatshap.vcf.gz
sample_whatshap.vcf.gz.tbi
```

VCF files should contain:

- phased genotypes (`0|1`, `1|0`)
- phase set (`PS`) information

Example directory structure:

```text
vcfs_wider/
├── sample1_whatshap.vcf.gz
├── sample1_whatshap.vcf.gz.tbi
├── sample2_whatshap.vcf.gz
└── sample2_whatshap.vcf.gz.tbi
```

---

## 2. Heterozygous carrier samples

Three heterozygous individuals carrying the variant or locus of interest must be provided.

These individuals are used to:

- identify shared cis-associated variants

- define the shared haplotypic interval

- reconstruct the candidate founder haplotype

The workflow searches for variants shared in at least two of the three heterozygous individuals to define the haplotype block boundaries.

Samples are provided as a comma-separated list:

```bash

-s 23GM4915,25_3624,AT015

```

Sample names must match the VCF filenames:

```text

23GM4915_whatshap.vcf.gz

25_3624_whatshap.vcf.gz

AT015_whatshap.vcf.gz

```

---

## 3. Genomic position of interest

Genomic coordinate associated with the target locus.

Example:

```bash
chr16:89511445
```

The workflow is compatible with:

- SNVs
- indels
- STR expansions
- structural variant-associated loci

---

## 4. Variant annotation file

Tabular annotation file containing genomic positions and rsIDs. You can download this file from the UCSC Table Browser [https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=varRep&hgta_track=dbSnp155Composite&hgta_table=dbSnp155].

It will be necessary to select a region of interest, be sure it is wide enough. In output format, "selected fields from primary and related tables" can be selected so only the next columns are exported: 
- 	chrom:	Reference sequence chromosome or scaffold
-	chromStart:	Start position in chrom
-	chromEnd:	End position in chrom
-	name:	dbSNP Reference SNP (rs) identifier
-	ref:	Reference allele; usually major allele, but may be minor allele
- 	alts:	Alternate (non-reference) alleles; may include major allele

Be sure it is tsv format. 

This file is used to annotate variants and identify known SNPs within the reconstructed haplotypes.

---

# Command-line arguments

The pipeline is executed using:

```bash
bash run_haplotyping.sh \
    -d <vcf_directory> \
    -p <genomic_position> \
    -s <sample1,sample2,...> \
    -v <variant_annotation_file> \
    -o <output_directory>
```

## Arguments

| Argument | Description |
|---|---|
| `-d` | Directory containing phased VCF files |
| `-p` | Genomic position of interest |
| `-s` | Comma-separated sample names |
| `-v` | Variant annotation file |
| `-o` | Output directory |

---

# Example execution

```bash
bash run_haplotyping.sh \
    -d vcfs \
    -p chr16:89511445 \
    -s 23GM4915,25_3624,AT015 \
    -v hg38_891-092M.summary.tsv \
    -o haplotyping
```

---

# Workflow

## Step 1 — Extract phase block-associated heterozygous variants

For each sample:

- identify the phase set (PS) linked to the locus of interest
- extract phased heterozygous variants belonging to the same haplotype block

Output:

```text
$output/txtfiles_heteros/
```

---

## Step 2 — Define shared haplotypic interval

Script:

```bash
python3 1-define_block.py
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
$output/txtfiles_all/
```

---

## Step 4 — Compare haplotypes across samples

Script:

```bash
python3 2-search_variants.py
```

This step identifies:

- shared phased variants
- common cis-associated alleles
- candidate shared haplotypes
- potential founder backgrounds

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
│   ├── 1-define_block.py
│   └── 2-search_variants.py
│
├── example/
│
└── README.md
```


---

# Citation

If you use this repository, please cite our work. 

