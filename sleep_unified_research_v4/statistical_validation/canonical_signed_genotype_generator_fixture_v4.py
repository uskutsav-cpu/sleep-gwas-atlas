#!/usr/bin/env python3
"""Tiny invented-genotype algebra fixture; cannot read assets or launch estimators."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).with_name("canonical_signed_genotype_generator_fixture_receipt_v4.json")


def matrix(correlations, h2):
    d = np.sqrt(np.asarray(h2, dtype=np.float64))
    return np.asarray(correlations, dtype=np.float64) * d[:, None] * d[None, :]


def main():
    if OUT.exists():
        raise SystemExit("preserve existing receipt; choose a distinct version")
    # Invented dosage matrix, samples by SNPs. No biological inputs are opened.
    dosage = np.array([[((s + 1) * (v + 3) + s * s + v * v) % 3
                        for v in range(12)] for s in range(7)], dtype=float)
    centered = dosage - dosage.mean(axis=0)
    keep = centered.std(axis=0) > 0
    G = centered[:, keep] / centered[:, keep].std(axis=0)
    n, m = G.shape
    R = G.T @ G / n
    scenarios = {
        "C01": (matrix([[1, .3, 0, 0], [.3, 1, 0, 0], [0, 0, 1, .3], [0, 0, .3, 1]], [.2, .1, .2, .1]), np.eye(4)),
        "C02": (matrix([[1, .3, .3], [.3, 1, .8], [.3, .8, 1]], [.1, .2, .2]), np.array([[1, 0, 0], [0, 1, .3], [0, .3, 1]])),
        "C03": (matrix([[1, .3, .3], [.3, 1, .8], [.3, .8, 1]], [.1, .2, .2]), np.array([[1, 0, 0], [0, 1, .3], [0, .3, 1]])),
        "C04": (matrix([[1, .5, .3], [.5, 1, .3], [.3, .3, 1]], [.1, .03, .1]), np.array([[1, .2, 0], [.2, 1, 0], [0, 0, 1]])),
    }
    checks = 0
    records = {}

    def close(a, b):
        nonlocal checks
        np.testing.assert_allclose(a, b, atol=1e-12, rtol=1e-12)
        checks += 1

    for name, (sigma_g, sigma_e) in scenarios.items():
        q = len(sigma_g)
        Lg, Le = np.linalg.cholesky(sigma_g), np.linalg.cholesky(sigma_e)
        close(Lg @ Lg.T, sigma_g)
        close(Le @ Le.T, sigma_e)
        # Fixed deterministic latent vectors test arithmetic, never empirical calibration.
        U = np.arange(m * q, dtype=float).reshape(m, q) / 23 - .7
        H = np.arange(n * q, dtype=float).reshape(n, q) / 29 - .4
        beta = U @ Lg.T / np.sqrt(m)
        eta = H @ Le.T
        direct_signal, direct_noise = R @ beta, G.T @ eta / np.sqrt(n)
        latent_genetic = sum(G[:, lo:lo + 3] @ beta[lo:lo + 3]
                             for lo in range(0, m, 3)) / n
        streamed_signal = np.concatenate([G[:, lo:lo + 3].T @ latent_genetic
                                           for lo in range(0, m, 3)])
        streamed_noise = np.concatenate([G[:, lo:lo + 3].T @ eta / np.sqrt(n)
                                          for lo in range(0, m, 3)])
        close(direct_signal, streamed_signal)
        close(direct_noise, streamed_noise)
        N = np.full((m, q), 250000.0)
        if name == "C03":
            N[:] = np.array([125000, 250000, 500000])[np.arange(m) % 3, None]
        if name == "C04":
            N[:, 1] = 85000
        close(np.sqrt(N) * direct_signal + direct_noise,
              np.sqrt(N) * streamed_signal + streamed_noise)
        # Analytic covariance of vec(eps) and vec(signal), SNP-major flattening.
        Ae = np.kron(G.T / np.sqrt(n), Le)
        Ag = np.kron(R / np.sqrt(m), Lg)
        close(Ae @ Ae.T, np.kron(R, sigma_e))
        close(Ag @ Ag.T, np.kron(R @ R, sigma_g / m))
        close(np.diag(R), np.ones(m))
        # An allele flip changes both genotype column and source-oriented beta row.
        flip = np.ones(m); flip[::3] = -1
        flipped_G = G * flip
        flipped_beta = beta * flip[:, None]
        close(flipped_G.T @ (flipped_G @ flipped_beta / n), direct_signal * flip[:, None])
        close(flipped_G.T @ eta / np.sqrt(n), direct_noise * flip[:, None])
        records[name] = {"genetic_covariance": sigma_g.tolist(), "noise_correlation": sigma_e.tolist(),
                         "genetic_min_eigenvalue": float(np.linalg.eigvalsh(sigma_g).min()),
                         "noise_min_eigenvalue": float(np.linalg.eigvalsh(sigma_e).min()),
                         "PSD_by_cholesky_without_smoothing": True}
    bad = np.array([[1., 2.], [2., 1.]])
    try:
        np.linalg.cholesky(bad)
    except np.linalg.LinAlgError:
        checks += 1
    else:
        raise AssertionError("indefinite model silently accepted")
    receipt = {"scope": "invented 7-sample genotype algebra only; no asset decode, random replicates, native imports, fits, covariance outcomes or P values",
               "checks": checks, "fixture_shape": list(G.shape), "numpy_version": np.__version__,
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "scenarios": records, "full_realistic_calibration_pass": False,
               "population_LD_or_reference_matching_established": False}
    OUT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"checks": checks, "scope": receipt["scope"], "full_realistic_calibration_pass": False}))


if __name__ == "__main__":
    main()
