import pandas as pd
import matplotlib.pyplot as plt
import sys


def plot_haplotype_blocks(df, pos_mark, region_start, region_end, outdir):	
	plt.figure(figsize=(10,3))

	plt.axvspan(region_start, region_end, alpha=0.2, color="red")

	for i, row in df.iterrows():

		plt.hlines(
			y=i,
			xmin=row["start"],
			xmax=row["end"],
			linewidth=4
		)

		plt.text(
			row["end"] + 5000,
			i,
			row["sample"],
			va='center'
		)

	plt.axvline(x=pos_mark, linestyle='--')

	plt.yticks(range(len(df)), df["sample"])

	plt.xlabel("Genomic position")
	plt.ylabel("Sample")
	plt.title("Overlap of haplotype blocks")

	plt.tight_layout()
	plt.savefig(f"{outdir}/block_plot.png", dpi=300)

if __name__ == "__main__":
	df = pd.read_csv(sys.argv[1], sep="\t")
	pos_mark = int(sys.argv[2].split(":")[1])
	region_start = int(sys.argv[3])
	region_end = int(sys.argv[4])
	plot_haplotype_blocks(df, pos_mark, region_start, region_end, sys.argv[5])