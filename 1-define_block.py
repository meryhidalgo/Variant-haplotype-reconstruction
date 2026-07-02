import pandas as pd
import json, sys, os


# esta función recoge un txt con las variantes de un individuo y devuelve un dataframe con la variante de interés 
# y las variantes cis a ella, es decir, las que están en el mismo haplotipo
def variantes_cis(archivo, var_interes_pos):
    df = pd.read_csv(archivo, sep="\t", header=None,
                     names=["chr","pos", "REF", "ALT", "GT","PS"])

    var_int = df[df["pos"] == var_interes_pos]
    if var_int.empty:
        print(f"Variante de interés no encontrada en {archivo}")
        return pd.DataFrame(), archivo  # Variante no encontrada

    # Sacar haplotipo y PS de la variante de interés
    gt = var_int["GT"].iloc[0].replace("/", "|").split("|")
    #ps = var_int["PS"].iloc[0] # ya he filtrado por PS en el awk, así que no es necesario
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
    var_interest = sys.argv[1]
    var_interest_chr = var_interest.split(":")[0]
    var_interest_pos = int(var_interest.split(":")[1])

    txt_dir = sys.argv[2]  # Directorio donde se encuentran los archivos txt
    txt_files = [os.path.join(txt_dir, f) for f in os.listdir(txt_dir) if f.endswith("_heteros.txt")]
    if len(txt_files) == 0:
        print(f"No se encontraron archivos txt en el directorio {txt_dir}")
        sys.exit(1)
    elif len(txt_files) < 3 or len(txt_files) > 3:
        print(f"Se encontró un número diferente a 3 archivos txt en el directorio {txt_dir}. Asegúrate de que solo haya los archivos de interés.")
        sys.exit(1)
    
    hap1 = variantes_cis(txt_files[0], var_interest_pos)
    hap2 = variantes_cis(txt_files[1], var_interest_pos)
    hap3 = variantes_cis(txt_files[2], var_interest_pos)

    dbSNP_file = sys.argv[3]
    dbSNP_all = pd.read_csv(dbSNP_file, sep="\t", header=0)

    hap1["sample"] = "hetero1"
    hap2["sample"] = "hetero2"
    hap3["sample"] = "hetero3"

    # Concatenar
    all_df = pd.concat([hap1, hap2, hap3])

    # Contar ocurrencias por variante
    counts = (
        all_df
        .groupby(["pos", "variant", "REF", "ALT"])["sample"]
        .nunique()
        .reset_index(name="n_samples")
    )

    shared = counts[counts["n_samples"] >= 2]

    # Merging with dbSNP_all to get additional information
    merged = pd.merge(shared, dbSNP_all, left_on="pos", right_on="end", how="left")

    found_inRS = pd.DataFrame(columns=["pos", "variant", "dbSNP_id", "ref", "alt"])
    for idx, row in merged.iterrows():
        alternatives = [a.strip() for a in str(row["alt"]).split(",")]
        if row["REF"] == row["ref"] and row["ALT"] in alternatives:
            found_inRS.loc[len(found_inRS)] = {"pos": row["pos"], "variant": row["variant"], "dbSNP_id": row["id"], "ref": row["ref"], "alt": row["alt"]}
        #else:
            #print(f"Posición {row['pos']} no coincide con dbSNP: REF {row['REF']} vs {row['ref']}, ALT {row['ALT']} vs {alternatives}")
    print(f"Se han encontrado {len(found_inRS)} variantes en dbSNP", file=sys.stderr)

    outdir = sys.argv[4]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    found_inRS.to_csv(
        f"{outdir}/heteros_cis_dbSNP_{len(found_inRS)}variants.txt",
        sep="\t",
        index=False
    )
    print(
        f"La región de interés está entre las posiciones "
        f"{found_inRS.iloc[0]['pos']} y {found_inRS.iloc[-1]['pos']}",
        file=sys.stderr
    )

    result = {
        "start": int(found_inRS.iloc[0]["pos"]),
        "end": int(found_inRS.iloc[-1]["pos"]),
        "n_variants": len(found_inRS),
        "hetero_cis_file": f"{outdir}/heteros_cis_dbSNP_{len(found_inRS)}variants.txt"
    }
    print(json.dumps(result))