#!/usr/bin/env python3
"""Generate kernel PCA reference fixture for R cross-validation.

Run from the Python conga repo root:
    python3 /path/to/rconga/tests/generate_kernel_pca_fixture.py

Reads:   tests/testthat/fixtures/dash.csv  (first 100 mouse TCRs)
Writes:  tests/testthat/fixtures/kernel_pca_ref.json

Uses sklearn KernelPCA(kernel='precomputed') as the reference
implementation, matching the logic in conga/preprocess.py lines 2203-2256.
"""
import csv
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, '.')

from sklearn.decomposition import KernelPCA
from conga.tcrdist.tcr_distances import TcrDistCalculator

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(SCRIPT_DIR, 'testthat', 'fixtures')
DASH_CSV    = os.path.join(FIXTURE_DIR, 'dash.csv')
OUTPUT_JSON = os.path.join(FIXTURE_DIR, 'kernel_pca_ref.json')

N_TCRS = 100
N_COMPONENTS = 20

# ---------------------------------------------------------------------------
# Load first 100 TCRs from DASH dataset
# ---------------------------------------------------------------------------
with open(DASH_CSV, 'r') as f:
    reader = csv.DictReader(f)
    rows = [next(reader) for _ in range(N_TCRS)]

print(f"Loaded {len(rows)} TCRs from {DASH_CSV}")

# Build TCR tuples in CoNGA format: ((va, ja, cdr3a), (vb, jb, cdr3b))
tcrs = []
for r in rows:
    tcr = (
        (r['v_a_gene'], r['j_a_gene'], r['cdr3_a_aa']),
        (r['v_b_gene'], r['j_b_gene'], r['cdr3_b_aa']),
    )
    tcrs.append(tcr)

n = len(tcrs)
assert n == N_TCRS

# ---------------------------------------------------------------------------
# Compute NxN TCRdist distance matrix
# ---------------------------------------------------------------------------
print(f"Computing {n}x{n} TCRdist matrix (mouse)...")
tdist = TcrDistCalculator('mouse')

start = time.time()
D = np.array([tdist(x, y) for x in tcrs for y in tcrs]).reshape((n, n))
elapsed = time.time() - start
print(f"Distance matrix computed in {elapsed:.1f}s")
print(f"D.max() = {D.max():.2f}, D.min() = {D.min():.2f}")

# ---------------------------------------------------------------------------
# Default kernel: gram = max(0, 1 - D / Dmax)
# ---------------------------------------------------------------------------
print(f"\nDefault kernel: n_components={N_COMPONENTS}")
Dmax = float(D.max())
gram_default = np.maximum(0.0, 1 - (D / Dmax))

pca_default = KernelPCA(kernel='precomputed', n_components=N_COMPONENTS)
xy_default = pca_default.fit_transform(gram_default)

eigenvalues_default = pca_default.eigenvalues_.tolist()
embeddings_default = xy_default.tolist()
total_variance_default = float(np.sum(pca_default.eigenvalues_))

print(f"  Dmax = {Dmax:.4f}")
print(f"  eigenvalues (first 5): {eigenvalues_default[:5]}")
print(f"  total_variance = {total_variance_default:.6f}")
print(f"  embeddings shape: {xy_default.shape}")

# ---------------------------------------------------------------------------
# Gaussian kernel: gram = exp(-0.5 * (D / sdev)^2)
# ---------------------------------------------------------------------------
GAUSSIAN_SDEV = 100.0
print(f"\nGaussian kernel: sdev={GAUSSIAN_SDEV}, n_components={N_COMPONENTS}")
gram_gaussian = np.exp(-0.5 * (D / GAUSSIAN_SDEV) ** 2)

pca_gaussian = KernelPCA(kernel='precomputed', n_components=N_COMPONENTS)
xy_gaussian = pca_gaussian.fit_transform(gram_gaussian)

eigenvalues_gaussian = pca_gaussian.eigenvalues_.tolist()
embeddings_gaussian = xy_gaussian.tolist()
total_variance_gaussian = float(np.sum(pca_gaussian.eigenvalues_))

print(f"  eigenvalues (first 5): {eigenvalues_gaussian[:5]}")
print(f"  total_variance = {total_variance_gaussian:.6f}")
print(f"  embeddings shape: {xy_gaussian.shape}")

# ---------------------------------------------------------------------------
# Write fixture JSON
# ---------------------------------------------------------------------------
fixture = {
    "n": n,
    "organism": "mouse",
    "n_components": N_COMPONENTS,
    "default_kernel": {
        "Dmax": Dmax,
        "eigenvalues": eigenvalues_default,
        "embeddings": embeddings_default,
        "total_variance": total_variance_default,
    },
    "gaussian_kernel": {
        "sdev": GAUSSIAN_SDEV,
        "eigenvalues": eigenvalues_gaussian,
        "embeddings": embeddings_gaussian,
        "total_variance": total_variance_gaussian,
    },
}

with open(OUTPUT_JSON, 'w') as f:
    json.dump(fixture, f, indent=2)

print(f"\nFixture written to {OUTPUT_JSON}")
print(f"  n: {n}")
print(f"  n_components: {N_COMPONENTS}")
print(f"  default eigenvalues: {len(eigenvalues_default)}")
print(f"  gaussian eigenvalues: {len(eigenvalues_gaussian)}")
