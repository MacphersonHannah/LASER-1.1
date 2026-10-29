#!/usr/bin/env python3

import argparse
import subprocess
import glob
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--barcode_list", required=True)
    parser.add_argument("--abundance_script", required=True)
    return parser.parse_args()

def run_abundance(prefix, cluster_report, abund_file, log_file, mapped_output):
    cmd = f"python {args.abundance_script} {prefix} {cluster_report}"
    subprocess.run(f"{cmd} > {abund_file} 2> {log_file}", shell=True, check=True)

    with open(abund_file) as fin, open(mapped_output, "w") as fout:
        for line in fin:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            fields = line.split()
            if len(fields) >= 3:
                fout.write(f"{fields[0]}\t{fields[1]}\n")

if __name__ == "__main__":
    args = parse_args()

    with open(args.barcode_list) as f:
        barcodes = f.read().strip().splitlines()

    for folder in barcodes + ["bulk"]:
        prefix = f"Collapsed_isoforms/{folder}/{folder}.collapsed"
        cluster_glob = glob.glob(f"Cluster/{folder}/*_clustered.cluster_report.csv")
        if not cluster_glob:
            raise FileNotFoundError(f"No cluster report for {folder}")
        cluster_report = cluster_glob[0]

        abund_file = f"Collapsed_isoforms/{folder}/{folder}.collapsed.abundance.txt"
        log_file = f"Collapsed_isoforms/{folder}/{folder}_get_abundance.log"
        mapped_output = f"Collapsed_isoforms/{folder}/{folder}_mapped_fl_count.txt"

        run_abundance(prefix, cluster_report, abund_file, log_file, mapped_output)
