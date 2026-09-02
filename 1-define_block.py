import pandas as pd
import json, sys, os


# This function retrieves a text file containing the variants of an individual and returns a dataframe with the variant of interest.
# And the cis variants of it, that is, those that are in the same haplotype.
def variantes_cis(archivo, var_interes_pos):
    df = pd.read_csv(archivo, sep="\t", header=None,
                     names=["chr","pos", "REF", "ALT", "GT","PS"])

    var_int = df[df["pos"] == var_interes_pos]
    if var_int.empty:
        print(f"Variant of interest not found in {archivo}")
        return pd.DataFrame(), archivo  # Variant not found

    # Extract haplotype and PS of the variant of interest
    gt = var_int["GT"].iloc[0].replace("/", "|").split("|")
    #ps = var_int["PS"].iloc[0] # I've already filtered by PS in awk, so it's not necessary.
    hap_int = 0 if gt[0] == "1" else 1
    #print(hap_int)
    results = pd.DataFrame(columns=["pos", "REF", "ALT", "variant"])
    for idx, row in df.iterrows():
        gt = row["GT"].replace("/", "|").split("|")
        if (gt[0] == gt[1] == 1):
            variant = row["ALT"]
        else:
            variant = row["REF"] if gt[hap_int] == "0" else row["ALT"]
        results.loc[len(results)] = {"pos": row["pos"], "REF": row["REF"], "ALT": row["ALT"], "variant": variant}
    return results



if __name__ == "__main__":
    min_samples = int(sys.argv[1])
    var_interest = sys.argv[2]
    var_interest_chr = var_interest.split(":")[0]
    var_interest_pos = int(var_interest.split(":")[1])

    txt_dir = sys.argv[3]  # Directory where the txt files are located
    txt_files = [os.path.join(txt_dir, f) for f in os.listdir(txt_dir) if f.endswith("_heteros.txt")]
    if len(txt_files) == 0:
        print(f"No txt files were found in the directory {txt_dir}", file=sys.stderr)
        sys.exit(1)
    elif len(txt_files) < min_samples:
        print(f"A number less than {min_samples} txt files was found in the directory {txt_dir}. Please adjust the number of files.", file=sys.stderr)
        sys.exit(1)
    

    dbSNP_file = sys.argv[4]
    if not os.path.exists(dbSNP_file):
        print(f"dbSNP file {dbSNP_file} does not exist.", file=sys.stderr)
        sys.exit(1)
    dbSNP_all = pd.read_csv(dbSNP_file, sep="\t", header=0)

    for f in txt_files:
        hap = variantes_cis(f, var_interest_pos)
        hap["sample"] = os.path.basename(f).replace("_heteros.txt", "")
        all_df = hap if 'all_df' not in locals() else pd.concat([all_df, hap])

    # Count occurrences by variant
    counts = (
        all_df
        .groupby(["pos", "variant", "REF", "ALT"])["sample"]
        .nunique()
        .reset_index(name="n_samples")
    )

    shared = counts[counts["n_samples"] >= min_samples]

    # Merging with dbSNP_all to get additional information
    merged = pd.merge(shared, dbSNP_all, left_on="pos", right_on="end", how="left")

    found_inRS = pd.DataFrame(columns=["pos", "variant", "dbSNP_id", "ref", "alt"])
    for idx, row in merged.iterrows():
        alternatives = [a.strip() for a in str(row["alt"]).split(",")]
        if row["REF"] == row["ref"] and row["ALT"] in alternatives:
            found_inRS.loc[len(found_inRS)] = {"pos": row["pos"], "variant": row["variant"], "dbSNP_id": row["id"], "ref": row["ref"], "alt": row["alt"]}
        #else:
            #print(f"Position {row['pos']} does not match dbSNP: REF {row['REF']} vs {row['ref']}, ALT {row['ALT']} vs {alternatives}")
    print(f"{len(found_inRS)} variants found in dbSNP", file=sys.stderr)

    outdir = sys.argv[5]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    found_inRS.to_csv(
        f"{outdir}/heteros_cis_dbSNP_{len(found_inRS)}variants.txt",
        sep="\t",
        index=False
    )
    print(
        f"The region of interest is between positions "
        f"{found_inRS.iloc[0]['pos']} and {found_inRS.iloc[-1]['pos']}",
        file=sys.stderr
    )

    result = {
        "start": int(found_inRS.iloc[0]["pos"]),
        "end": int(found_inRS.iloc[-1]["pos"]),
        "n_variants": len(found_inRS),
        "hetero_cis_file": f"{outdir}/heteros_cis_dbSNP_{len(found_inRS)}variants.txt"
    }
    print(json.dumps(result))