import pandas as pd
import sys, os

def extract_ps_from_txtfile(path, pos):
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "pos", "ref", "alt", "gt", "ps"])
    hit = df[df["pos"] == pos]
    if hit.empty:
        return None
    ps_values = hit["ps"].unique()
    # si hay más de uno o missing
    if len(ps_values) == 0 or (len(ps_values) == 1 and ps_values[0] == "."):
        return None
    return ps_values[0]

def extract_phased_alleles(variant_file, block_id):
    cols = ["chrom", "pos", "ref", "alt", "gt", "phase"]
    df = pd.read_csv(variant_file, sep="\t", header=None, names=cols)
    records = []
    for _, row in df.iterrows():
        info = "No same block"
        g1, g2 = row["gt"].split("|")
        hap1 = row["ref"] if g1 == "0" else row["alt"]
        hap2 = row["ref"] if g2 == "0" else row["alt"]
        if row["gt"] == "1|1":
            info = "homo"
        elif block_id != "." and str(row["phase"]) == block_id:
            info = "Same block"
        records.append({
            "pos": row["pos"],
            "hap1": hap1,
            "hap2": hap2,
            "info": info
        })
    result = pd.DataFrame(records)
    return result


if __name__ == "__main__":
    var_interest = sys.argv[1]
    var_interest_chr = var_interest.split(":")[0]
    var_interest_pos = int(var_interest.split(":")[1])

    found_inRS_file = sys.argv[2]
    found_inRS= pd.read_csv(found_inRS_file, sep="\t")

    txt_dir = sys.argv[3]  # Directorio donde se encuentran los archivos txt
    txt_files = [os.path.join(txt_dir, f) for f in os.listdir(txt_dir) if f.endswith(".txt")]

    phased_all = {}

    for f in txt_files:
        ps = extract_ps_from_txtfile(f, var_interest_pos)
        sample = os.path.basename(f).replace(".txt", "")
        if ps is None:
            print(f"[WARN] {sample}: PS missing at pos {var_interest}")
            ps = "."
        #print(f"{sample}\t{var_interest}\t{ps}")
        phased_all[sample] = extract_phased_alleles(f, ps)


    merged = found_inRS.copy()
    for sample,df in phased_all.items():
        merged = pd.merge(df, merged, on="pos", how="right", suffixes=(f"_{sample}", ""))
        merged.rename(columns={f"hap1": f"hap1_{sample}", f"hap2": f"hap2_{sample}", f"info": f"info_{sample}"}, inplace=True)
        cols = [f"hap1_{sample}", f"hap2_{sample}"]
        mask = merged[cols].isna().any(axis=1)
        merged.loc[mask, cols] = merged.loc[mask, cols].apply(
            lambda col: col.fillna(merged.loc[mask, "ref"])
        )

        merged.loc[mask, f"info_{sample}"] = "not found"

        #merged = merged.apply(lambda row: row.fillna(row["ref"]) and row["info"]=="not found", axis=1)

    first = ["dbSNP_id", "pos", 'variant', 'ref', 'alt']
    merged = merged[first + [c for c in merged.columns if c not in first]]

    outdir = sys.argv[4]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    merged.to_csv(f"{outdir}/phased_allsamples_output_info.tsv", index=False, sep="\t")
    
    print(f"\nDONE! Output saved to {outdir}/phased_allsamples_output_info.tsv")