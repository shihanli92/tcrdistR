#!/usr/bin/env python3
"""
Benchmark Python tcrdist3 for comparison with tcrdistR.

Usage: python benchmarks/benchmark_vs_tcrdist3.py [path_to_dash.csv]

Requires: pip install tcrdist3 pandas
"""

import sys
import time
import pandas as pd
from tcrdist.repertoire import TCRrep


def main():
    if len(sys.argv) >= 2:
        dash_file = sys.argv[1]
    else:
        dash_file = "tests/testthat/fixtures/dash.csv"

    print(f"Loading data from: {dash_file}")
    dash = pd.read_csv(dash_file)
    dash = dash.rename(columns={
        "v_a_gene": "va",
        "cdr3_a_aa": "cdr3a",
        "v_b_gene": "vb",
        "cdr3_b_aa": "cdr3b",
    })
    print(f"Total TCRs: {len(dash)}\n")

    sizes = [50, 100, 200, 500, 1000, len(dash)]
    sizes = [s for s in sizes if s <= len(dash)]

    results = []
    for n in sizes:
        sub = dash.head(n).copy().reset_index(drop=True)
        npairs = n * (n - 1) / 2

        # tcrdist3 computes the matrix via TCRrep
        sub_rep = sub[["va", "cdr3a", "vb", "cdr3b"]].copy()
        sub_rep["count"] = 1

        t0 = time.time()
        tr = TCRrep(
            cell_df=sub_rep,
            organism="mouse",
            chains=["alpha", "beta"],
            compute_distances=True,
        )
        elapsed = time.time() - t0

        rate = npairs / elapsed if elapsed > 0 else float("inf")
        results.append((n, npairs, elapsed, rate))
        print(f"  n={n:5d} | {npairs:10.0f} pairs | {elapsed:7.4f} sec | {rate:10.0f} pairs/sec")

    print("\n--- tcrdist3 (Python) Benchmark Results ---")
    print(f"{'n':>6}  {'pairs':>12}  {'time_sec':>10}  {'pairs/sec':>12}")
    for n, npairs, elapsed, rate in results:
        print(f"{n:6d}  {npairs:12.0f}  {elapsed:10.4f}  {rate:12.0f}")

    import tcrdist
    print(f"\nPython version: {sys.version.split()[0]}")
    print(f"tcrdist3 version: {tcrdist.__version__}")


if __name__ == "__main__":
    main()
