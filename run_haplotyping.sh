#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage:"
    echo "  bash run_haplotyping.sh -d <vcf_directory> -p <position> -s <sample1,sample2,...> -v <variant-annotation-file> -o <output_directory>"
    exit 1
}

while getopts "d:p:s:v:o:" opt; do
    case $opt in
        d) VCF_DIR="$OPTARG" ;;
        p) POS="$OPTARG" ;;
        s) SAMPLE_STRING="$OPTARG" ;;
        v) VARIANT_ANNOTATION_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "${VCF_DIR:-}" || -z "${POS:-}" || -z "${SAMPLE_STRING:-}" ]]; then
    usage
fi

IFS=',' read -ra samples <<< "$SAMPLE_STRING"

mkdir -p ${OUTPUT_DIR}/txtfiles_heteros
mkdir -p ${OUTPUT_DIR}/txtfiles_all

echo "VCF directory: $VCF_DIR"
echo "Position: $POS"
echo "Samples: ${samples[*]}"
echo "Variant annotation file: $VARIANT_ANNOTATION_FILE"
echo "Output directory: $OUTPUT_DIR"

# bgzip and index VCFs if needed
for vcf in ${VCF_DIR}/*.vcf; do
    if [[ ! -f "${vcf}.gz" ]]; then
        bgzip -c "$vcf" > "${vcf}.gz"
        tabix -p vcf "${vcf}.gz"
    fi
done

# summary file for haplotype blocks
BLOCK_SUMMARY="${OUTPUT_DIR}/haplotype_blocks.tsv"

echo -e "sample\tBlock_id\tstart\tend" > "$BLOCK_SUMMARY"

# extract heterozygous phased variants
for sample in "${samples[@]}"; do
    VCF="${VCF_DIR}/${sample}_whatshap.vcf.gz"
    if [[ ! -f "$VCF" ]]; then
        echo "ERROR: VCF not found for sample $sample"
        continue
    fi
    echo "Processing $sample"
    PS=$(bcftools query -r "$POS" -f '[%PS\n]' "$VCF")
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' \
        "$VCF" |
    awk -v ps="$PS" '
        ($6 == ps && $5 ~ /\|/) {
            gsub(/\//,"|",$5);
            print
        }
    ' > "${OUTPUT_DIR}/txtfiles_heteros/${sample}_heteros.txt"

    # define block coordinates
    s=$(awk 'NR==1 {print $2}' "${OUTPUT_DIR}/txtfiles_heteros/${sample}_heteros.txt")
    e=$(awk '
        BEGIN{max=0}
        {
            if ($2 > max)
                max=$2
        }
        END{print max}
    ' "${OUTPUT_DIR}/txtfiles_heteros/${sample}_heteros.txt")
    echo -e "${sample}\t${PS}\t${s}\t${e}" >> "$BLOCK_SUMMARY"
done


# define shared block
json=$(python 1-define_block.py \
    "$POS" \
    "${OUTPUT_DIR}/txtfiles_heteros" \
    "$VARIANT_ANNOTATION_FILE" \
    "${OUTPUT_DIR}/")

start=$(echo "$json" | jq -r .start)
end=$(echo "$json" | jq -r .end)
hetero_cis_file=$(echo "$json" | jq -r .hetero_cis_file)
echo "Shared block:"
echo "START=$start"
echo "END=$end"

mkdir -p "${OUTPUT_DIR}/plots"

python 2-plot_block.py \
    "${OUTPUT_DIR}/haplotype_blocks.tsv" \
    "$POS" \
    "$start" \
    "$end" \
    "${OUTPUT_DIR}/plots/"

# extract all variants within interval
for vcf in ${VCF_DIR}/*.vcf.gz; do

    vcf_name=$(basename "$vcf" _whatshap.vcf.gz)
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' \
        "$vcf" |
    awk -v start="$start" -v end="$end" '
        BEGIN {OFS="\t"}
        ($1 == "chr16" && $2 >= start && $2 <= end) {
            gsub(/\//,"|",$5);
            print
        }
    ' > "${OUTPUT_DIR}/txtfiles_all/${vcf_name}.txt"
done

# compare haplotypes
python 3-search_variants.py \
    "$POS" \
    "$hetero_cis_file" \
    "${OUTPUT_DIR}/txtfiles_all" \
    "${OUTPUT_DIR}"