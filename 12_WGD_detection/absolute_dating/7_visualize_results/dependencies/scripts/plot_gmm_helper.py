#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.mixture import GaussianMixture


MIN_AGE = 5
MAX_AGE = 100


def run_plot_gmm(ages_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for age_file in sorted(ages_dir.glob("*.pep")):
        with open(age_file, "r", encoding="utf-8") as handle:
            data = [float(line.strip()) for line in handle if line.strip() and MIN_AGE <= float(line.strip()) <= MAX_AGE]
        if not data:
            continue
        values = np.array(data).reshape(-1, 1)
        lowest_bic = np.inf
        best_gmm = None
        for component_count in range(1, 4):
            gmm = GaussianMixture(n_components=component_count)
            gmm.fit(values)
            bic = gmm.bic(values)
            if bic < lowest_bic:
                lowest_bic = bic
                best_gmm = gmm
        plt.figure(figsize=(5, 3))
        plt.hist(values, bins=30, density=True, alpha=0.6, color="#f5f5f5", edgecolor="black")
        x_values = np.linspace(min(data), max(data), 1000).reshape(-1, 1)
        for mean, covariance, weight in zip(
            best_gmm.means_.flatten(),
            np.sqrt(best_gmm.covariances_).flatten(),
            best_gmm.weights_.flatten(),
        ):
            curve = weight * (1 / (covariance * np.sqrt(2 * np.pi))) * np.exp(-((x_values - mean) ** 2) / (2 * covariance**2))
            plt.plot(x_values, curve, linewidth=2, color="#e4a19f")
        plt.tick_params(axis="both", which="major", labelsize=11)
        plt.savefig(output_dir / f"{age_file.stem}_gmm.png", transparent=True, dpi=300, bbox_inches="tight")
        plt.close()
