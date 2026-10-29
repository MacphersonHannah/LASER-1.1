#!/usr/bin/env python3
import sys
import pandas as pd

# Fixed category orders
structural_order = [
    "novel_not_in_catalog",
    "incomplete-splice_match",
    "full-splice_match",
    "novel_in_catalog",
    "genic",
    "fusion"
]

subcategory_order = [
    "at_least_one_novel_splicesite",
    "intron_retention",
    "5prime_fragment",
    "alternative_3end5end",
    "3prime_fragment",
    "mono-exon",
    "alternative_3end",
    "multi-exon",
    "combination_of_known_junctions",
    "mono-exon_by_intron_retention",
    "combination_of_known_splicesites",
    "internal_fragment",
    "alternative_5end",
    "reference_match"
]

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <gene_id> <classification.txt> [ORFannotate_summary.tsv]")
        sys.exit(1)

    gene_id = sys.argv[1]
    sqanti_file = sys.argv[2]
    orf_file = sys.argv[3] if len(sys.argv) > 3 else None

    # SQANTI processing  
    df = pd.read_csv(sqanti_file, sep="\t", dtype=str)
    total_rows = df.shape[0]

    filt = df[(df["strand"] == "-") & (df["associated_gene"].str.contains(gene_id, na=False))]
    total_isoforms = filt.shape[0]

    struct_counts = filt["structural_category"].value_counts().to_dict()
    subcat_counts = filt["subcategory"].value_counts().to_dict()

    # Output  
    print(f"Input_file\t{sqanti_file}\n")
    print(f"Total_rows_pre-filtering\t{total_rows}")
    print(f"Total_unique_{gene_id}_isoforms\t{total_isoforms}\n")

    print("Counts_by_Structural_Category")
    for cat in structural_order:
        print(f"{cat}\t{struct_counts.get(cat, 0)}")
    for cat in struct_counts:
        if cat not in structural_order:
            print(f"{cat}\t{struct_counts[cat]}")

    print("\nCounts_by_Subcategory")
    for sub in subcategory_order:
        print(f"{sub}\t{subcat_counts.get(sub, 0)}")
    for sub in subcat_counts:
        if sub not in subcategory_order:
            print(f"{sub}\t{subcat_counts[sub]}")

    # Optional ORFannotate processing  
    if orf_file:
        orf_df = pd.read_csv(orf_file, sep="\t", dtype=str)
        # Match on transcript_id in ORFannotate and isoform in SQANTI
        matching_ids = set(filt["isoform"])
        orf_filt = orf_df[orf_df["transcript_id"].isin(matching_ids)]
        
        total_orf_isoforms = orf_filt.shape[0]
        coding_counts = orf_filt["coding_class"].value_counts().to_dict()

        print(f"\nTotal_filtered_unique_transcripts_ORFannotate\t{total_orf_isoforms}")
        print("Counts_by_Coding_Class")
        for cls in ["coding", "noncoding"]:
            print(f"{cls}\t{coding_counts.get(cls, 0)}")
        for cls in coding_counts:
            if cls not in ["coding", "noncoding"]:
                print(f"{cls}\t{coding_counts[cls]}")

if __name__ == "__main__":
    main()
