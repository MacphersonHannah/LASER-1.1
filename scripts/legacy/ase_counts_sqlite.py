#!/usr/bin/env python3
import argparse
import sqlite3
import pandas as pd
from scipy.stats import binomtest

print("[DEBUG] something", flush=True)

def parse_args():
    ap = argparse.ArgumentParser(description="Summarise ASE counts per gene/sample/haplotype (with debug)")
    ap.add_argument("--sqlite", required=True)
    ap.add_argument("--passes", required=True)
    ap.add_argument("--gene_bed", required=True)
    ap.add_argument("--sample", required=False,
                    help="If provided, restrict to this sample; otherwise compute all samples in DB")
    ap.add_argument("--out", required=True)
    return ap.parse_args()


def main():
    args = parse_args()
    print("========== ASE DEBUG MODE ==========")
    print(f"[DEBUG] SQLite DB: {args.sqlite}")
    print(f"[DEBUG] PASS CSV: {args.passes}")
    print(f"[DEBUG] BED file: {args.gene_bed}")
    print(f"[DEBUG] Sample (if provided): {args.sample}")

    
    # 1. Read gene BED and build mapping tables
    
    map_df = pd.read_csv(
        args.gene_bed, sep="\t", header=None,
        names=["chrom", "start", "end", "gene_name", "gene_id"]
    )
    print(f"[DEBUG] Loaded BED with {len(map_df)} rows.")
    print("[DEBUG] BED example rows:")
    print(map_df.head(3).to_string(index=False))

    name_to_id = dict(zip(map_df["gene_name"], map_df["gene_id"]))
    id_to_name = dict(zip(map_df["gene_id"], map_df["gene_name"]))

    
    # 2. Read passes table (human-readable gene names)
    
    passes = pd.read_csv(args.passes)
    print(f"[DEBUG] Loaded PASS table with shape {passes.shape}")
    print(f"[DEBUG] PASS columns: {list(passes.columns)}")
    print(f"[DEBUG] Example PASS entries:\n{passes.head(3).to_string(index=False)}")

    pass_cols = [c for c in passes.columns if c != "gene"]

    # Translate gene names -> Ensembl IDs
    passes["gene_id"] = passes["gene"].map(name_to_id)
    unmapped = passes[passes["gene_id"].isna()]
    if not unmapped.empty:
        print(f"[DEBUG] WARNING: {len(unmapped)} PASS genes couldn't be mapped to Ensembl IDs.")
        print(f"[DEBUG] Example unmapped genes: {unmapped['gene'].head(5).tolist()}")

    all_gene_ids = passes["gene_id"].dropna().unique().tolist()
    print(f"[DEBUG] Total unique mapped gene IDs: {len(all_gene_ids)}")
    print(f"[DEBUG] Example gene IDs for query: {all_gene_ids[:5]}")

    
    # 3. Open SQLite database
    
    conn = sqlite3.connect(args.sqlite)
    cur = conn.cursor()

    # Detect which schema version is present
    tables = [r[0] for r in cur.execute("SELECT name FROM sqlite_master WHERE type='table'")]
    print(f"[DEBUG] Tables in SQLite: {tables}")

    if "reads" in tables and "transcripts" in tables:
        schema_type = "normalized"
        db_samples = [r[0] for r in cur.execute("SELECT DISTINCT sample FROM reads")]
    elif "isoform_lookup" in tables:
        schema_type = "legacy"
        db_samples = [r[0] for r in cur.execute("SELECT DISTINCT sample FROM isoform_lookup")]
    else:
        raise RuntimeError(f"[ERROR] No recognised tables found in {args.sqlite}")

    print(f"[DEBUG] Detected schema: {schema_type}")
    print(f"[DEBUG] Found {len(db_samples)} distinct samples in DB: {db_samples}")

    # 4. Fetch all rows for all PASS genes
    
    ph = ",".join("?" * len(all_gene_ids)) if all_gene_ids else "NULL"
    query = f"""
        SELECT r.sample, t.sqanti_associated_gene AS gene, r.haplotype
        FROM reads r
        JOIN transcripts t ON r.pbid = t.pbid
        WHERE t.sqanti_associated_gene IN ({ph})
        """

    print("[DEBUG] Executing SQL query for PASS genes...")
    df = (
        pd.read_sql_query(query, conn, params=all_gene_ids)
        if all_gene_ids else
        pd.DataFrame(columns=["sample", "gene", "haplotype"])
    )
    conn.close()
    print(f"[DEBUG] Retrieved {len(df)} rows from DB.")
    print(f"[DEBUG] Example SQL output:\n{df.head(5).to_string(index=False)}")

    if df.empty:
        open(args.out, "w").close()
        print("[WARN] No matching reads found in database for given genes.")
        print("[DEBUG] Possible reasons:")
        print("  - sqanti_associated_gene uses Ensembl IDs but your PASS table used names (or vice versa)")
        print("  - mismatch between BED and DB gene identifiers")
        return

    df["haplotype"] = df["haplotype"].fillna("unphased").astype(str)



    # 5. Process per sample
    
    # Decide which samples to process
    if args.sample:
        samples_to_process = [args.sample]
        print(f"[DEBUG] Restricting to single sample: {args.sample}")
    else:
        samples_to_process = db_samples
        print(f"[DEBUG] Processing all {len(samples_to_process)} samples from DB.")

    results = []
    for s in samples_to_process:
        print(f"[DEBUG] Processing sample: {s}")

        if s not in passes.columns:
            print(f"[WARN] Sample {s} not in PASSES columns; skipping.")
            continue

        pass_ids = passes.loc[
            passes[s].astype(str).str.upper() == "PASS", "gene_id"
        ].dropna().unique().tolist()
        print(f"[DEBUG] {s}: {len(pass_ids)} PASS genes mapped to Ensembl IDs.")
        if not pass_ids:
            continue

        sdf = df[(df["sample"] == s) & (df["gene"].isin(pass_ids))]
        print(f"[DEBUG] {s}: {len(sdf)} matching records in SQL results.")

        for gene, gdf in sdf.groupby("gene"):
            counts = gdf["haplotype"].value_counts().to_dict()
            h1, h2 = counts.get("1", 0), counts.get("2", 0)
            total = h1 + h2
            pval = binomtest(h1, total, 0.5).pvalue if total > 0 else None

            results.append({
                "sample": s,
                "gene_id": gene,
                "gene_name": id_to_name.get(gene, "NA"),
                "H1": h1,
                "H2": h2,
                "UNPHASED": counts.get("unphased", 0),
                "TOTAL": total,
                "BINOM_P": pval
            })

        print(f"[DEBUG] {s}: Added {len(results)} total rows so far.")

    
    # 6. Write output
    
    out_df = pd.DataFrame(results)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"[INFO] Wrote {len(out_df)} rows to {args.out}")
    print("========== DEBUG RUN COMPLETE ==========")


if __name__ == "__main__":
    main()
