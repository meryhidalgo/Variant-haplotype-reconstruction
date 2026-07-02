#!/bin/bash
#conda activate haplotyping 
# conda create -n haplotyping -c conda-forge -c bioconda bcftools htslib tabix pandas

#define the heteros of interest
samples=("23GM4915" "25_3624" "AT015")

#define the position of interest
POS="chr16:89511445"

#bgzipped vcfs if not
for vcf in vcfs_wider/*.vcf; do
	if [[ ! -f "${vcf}.gz" ]]; then
		bgzip -c $vcf > ${vcf}.gz
		tabix -p vcf ${vcf}.gz
	fi
done

VCF1=vcfs_wider/${samples[1]}_whatshap.vcf.gz
VCF2=vcfs_wider/${samples[2]}_whatshap.vcf.gz
VCF3=vcfs_wider/${samples[3]}_whatshap.vcf.gz

mkdir -p haplotyping/txtfiles_heteros

VCF1_name=$(basename $VCF1 _whatshap.vcf.gz)
VCF2_name=$(basename $VCF2 _whatshap.vcf.gz)
VCF3_name=$(basename $VCF3 _whatshap.vcf.gz)

# SPG7 in chr16:89511445, expansion in chr14:23321471
PS=$(bcftools query -r $POS -f '[%PS\n]' $VCF1) #esperando que siempre se reporte en la posición 23321471
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' $VCF1 | awk -v ps="$PS" '($6 == ps && $5 ~ /\|/) {gsub(/\//,"|",$5); print}' > haplotyping/txtfiles_heteros/${VCF1_name}_heteros.txt

PS=$(bcftools query -r $POS -f '[%PS\n]' $VCF2) #esperando que siempre se reporte en la posición 23321471
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' $VCF2 | awk -v ps="$PS" '($6 == ps && $5 ~ /\|/)  {gsub(/\//,"|",$5); print}' > haplotyping/txtfiles_heteros/${VCF2_name}_heteros.txt

PS=$(bcftools query -r $POS -f '[%PS\n]' $VCF3) #esperando que siempre se reporte en la posición 23321471
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' $VCF3 | awk -v ps="$PS" '($6 == ps && $5 ~ /\|/)  {gsub(/\//,"|",$5); print}' > haplotyping/txtfiles_heteros/${VCF3_name}_heteros.txt

#python3 5-1-define_block.py POSITION TXT_DIRECTORY dbSNP_file outdir

json=$(python3 5-1-define_block.py "$POS" "haplotyping/txtfiles_heteros" "rs_search/all/hg38_891-092M.summary.tsv" "haplotyping")
start=$(echo $json | jq -r .start)
end=$(echo $json | jq -r .end)
hetero_cis_file=$(echo $json | jq -r .hetero_cis_file)


# aquí hay que definir la región de interés, que es la que contiene las variantes para 2/3 heteros
#antes estaba 89017370-90057865 y ahora 89243866-90008504
mkdir -p haplotyping/txtfiles_all

for vcf in vcfs_wider/*.vcf.gz; do
	vcf_name=$(basename $vcf _whatshap.vcf.gz)
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' $vcf | awk -v start=$start -v end=$end 'BEGIN{OFS="\t"} ($1 == "chr16" && $2 >= start && $2 <= end) {gsub(/\//,"|",$5); print}' > haplotyping/txtfiles_all/${vcf_name}.txt
done

#python3 5-2-search_variants.py POS hetero_cis_file TXT_DIRECTORY dbSNP_file outdir
python3 5-2-search_variants.py "$POS" $hetero_cis_file "haplotyping/txtfiles_all" "haplotyping"
