# Benchmarks: tcrdistR vs Python tcrdist3

Comparison of pairwise TCRdist matrix computation times using the DASH dataset
(1,924 mouse TCRs).

## Running

### R (tcrdistR)

```bash
cd tcrdistR/
Rscript benchmarks/benchmark_vs_tcrdist3.R
```

### Python (tcrdist3)

```bash
pip install tcrdist3 pandas
cd tcrdistR/
python benchmarks/benchmark_vs_tcrdist3.py
```

Both scripts use `tests/testthat/fixtures/dash.csv` by default. Pass a custom
path as the first argument if needed.

## What is measured

Both scripts compute the full N x N pairwise TCRdist matrix for paired
alpha/beta TCRs at sizes N = 50, 100, 200, 500, 1000, and 1924. Wall-clock
time is reported.
