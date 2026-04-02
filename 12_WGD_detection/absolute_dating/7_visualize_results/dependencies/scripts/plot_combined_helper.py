#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from matplotlib.colors import to_rgb


INPUT_SUFFIX = ".pep"
MIN_AGE = 5
MAX_AGE = 80
BIN_COUNT = 48
MIN_POINTS_PER_SAMPLE = 5
CLADE_COLORS = {
    "Clade_5": "#c45a4a",
    "Clade_4": "#5b6fb8",
    "Clade_3": "#4b9bb0",
    "Clade_2": "#5da97a",
    "Clade_1": "#8a80c8",
    "Other": "#666666",
}
SMOOTH_KERNEL = np.array([1, 2, 3, 2, 1], dtype=float)
HEATMAP_GAMMA = 1.8
SAMPLE_PANEL_COLUMNS = 6


def read_sample_ages(file_path: Path) -> list[float]:
    with open(file_path, "r", encoding="utf-8") as handle:
        return [value for value in (float(line.strip()) for line in handle if line.strip()) if MIN_AGE <= value <= MAX_AGE]


def load_all_samples(ages_dir: Path):
    sample_data = []
    for each_file in sorted(ages_dir.glob(f"*{INPUT_SUFFIX}")):
        ages = read_sample_ages(each_file)
        if len(ages) >= MIN_POINTS_PER_SAMPLE:
            sample_data.append((each_file.name, np.array(ages, dtype=float)))
    return sample_data


def read_clade_groups(clade_file: Path):
    groups = []
    current = []
    for line in clade_file.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text:
            continue
        if text == "==":
            if current:
                groups.append(current)
                current = []
            continue
        current.append(text)
    if current:
        groups.append(current)
    return groups


def build_sample_to_clade(clade_groups):
    sample_to_clade = {}
    for index, group in enumerate(clade_groups):
        clade_name = f"Clade_{index + 1}"
        for sample_name in group:
            sample_to_clade[sample_name] = clade_name
    return sample_to_clade


def load_display_name_map(mapping_file: Path | None):
    if mapping_file is None or not mapping_file.exists():
        return {}
    result = {}
    with mapping_file.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            result[row["current_sample_name"]] = row["standard_latin_name"]
    return result


def read_tree_leaf_order(tree_file: Path):
    tree = Phylo.read(tree_file, "newick")
    tree.ladderize()
    return tree, [terminal.name for terminal in tree.get_terminals()]


def reorder_sample_data_by_tree(sample_data, leaf_order):
    sample_map = {sample_name: ages for sample_name, ages in sample_data}
    ordered = [(sample_name, sample_map[sample_name]) for sample_name in leaf_order if sample_name in sample_map]
    missing = [sample_name for sample_name in leaf_order if sample_name not in sample_map]
    return ordered, missing


def smooth_curve(values):
    kernel = SMOOTH_KERNEL / SMOOTH_KERNEL.sum()
    return np.convolve(values, kernel, mode="same")


def compute_clade_densities(sample_data, sample_to_clade, bin_edges):
    clade_ages = {}
    for sample_name, ages in sample_data:
        clade_name = sample_to_clade.get(sample_name, "Other")
        clade_ages.setdefault(clade_name, []).extend(ages.tolist())
    densities = {}
    for clade_name, ages in clade_ages.items():
        densities[clade_name] = np.histogram(np.array(ages), bins=bin_edges, density=True)[0]
    return densities


def build_heatmap(sample_data, bin_edges, sample_to_clade):
    rows, names, sizes, peaks, clades = [], [], [], [], []
    for sample_name, ages in sample_data:
        counts, _ = np.histogram(ages, bins=bin_edges)
        if counts.sum() == 0:
            continue
        normalized = counts / counts.sum()
        peak_index = int(np.argmax(normalized))
        rows.append(normalized)
        names.append(sample_name)
        sizes.append(len(ages))
        peaks.append((bin_edges[peak_index] + bin_edges[peak_index + 1]) / 2)
        clades.append(sample_to_clade.get(sample_name, "Other"))
    return np.array(rows), names, sizes, peaks, clades


def build_colored_heatmap(heatmap, ordered_clades):
    rgb_heatmap = np.ones((heatmap.shape[0], heatmap.shape[1], 3), dtype=float)
    for row_index, clade_name in enumerate(ordered_clades):
        base_color = np.array(to_rgb(CLADE_COLORS.get(clade_name, "#666666")))
        row_values = heatmap[row_index]
        row_low = np.percentile(row_values, 35)
        row_high = np.percentile(row_values, 95)
        row_values = np.clip((row_values - row_low) / (row_high - row_low), 0, 1) if row_high > row_low else np.zeros_like(row_values)
        row_values = np.clip((row_values - 0.1) / 0.9, 0, 1)
        enhanced = np.power(row_values, HEATMAP_GAMMA)
        rgb_heatmap[row_index, :, :] = 1 - (1 - base_color) * enhanced[:, None]
    return rgb_heatmap


def draw_tree_axis(ax, tree, leaf_order, display_name_map):
    leaf_y = {sample_name: index + 0.5 for index, sample_name in enumerate(leaf_order)}
    y_cache = {}
    depth_map = tree.depths() or tree.depths(unit_branch_lengths=True)
    max_depth = max(depth_map[terminal] for terminal in tree.get_terminals())

    def node_y(clade):
        if clade in y_cache:
            return y_cache[clade]
        y_value = leaf_y[clade.name] if clade.is_terminal() else float(sum(node_y(child) for child in clade.clades)) / len(clade.clades)
        y_cache[clade] = y_value
        return y_value

    def draw_clade(clade, parent_x):
        x_value = depth_map[clade]
        y_value = node_y(clade)
        if parent_x is not None:
            ax.plot([parent_x, x_value], [y_value, y_value], color="black", linewidth=0.7)
        if clade.clades:
            child_ys = [node_y(child) for child in clade.clades]
            ax.plot([x_value, x_value], [min(child_ys), max(child_ys)], color="black", linewidth=0.7)
            for child in clade.clades:
                draw_clade(child, x_value)
        else:
            ax.text(max_depth + 0.45, y_value, display_name_map.get(clade.name, clade.name), va="center", ha="left", fontsize=6.2, clip_on=False)

    draw_clade(tree.root, None)
    ax.set_xlim(0, max_depth + 4.0)
    ax.set_ylim(len(leaf_order), 0)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_sample_line_panels(sample_data, sample_to_clade, display_name_map, bin_edges, output_dir):
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    sample_curves = []
    max_curve_value = 0.0
    for sample_name, ages in sample_data:
        counts, _ = np.histogram(ages, bins=bin_edges, density=True)
        curve = smooth_curve(counts)
        sample_curves.append({"sample_name": sample_name, "clade_name": sample_to_clade.get(sample_name, "Other"), "curve": curve, "count": len(ages), "peak_index": int(np.argmax(curve))})
        max_curve_value = max(max_curve_value, float(np.max(curve)))
    if not sample_curves:
        return
    panel_count = len(sample_curves)
    ncols = min(SAMPLE_PANEL_COLUMNS, panel_count)
    nrows = math.ceil(panel_count / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 3.2, nrows * 2.3), sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    y_limit = max_curve_value * 1.08 if max_curve_value > 0 else 1.0
    for ax, curve_info in zip(axes, sample_curves):
        base_color = CLADE_COLORS.get(curve_info["clade_name"], "#666666")
        sample_name = curve_info["sample_name"]
        curve = curve_info["curve"]
        peak_x = bin_centers[curve_info["peak_index"]]
        peak_y = curve[curve_info["peak_index"]]
        ax.plot(bin_centers, curve, color=base_color, linewidth=1.8)
        ax.fill_between(bin_centers, curve, color=base_color, alpha=0.18)
        ax.scatter([peak_x], [peak_y], s=9, color=base_color, zorder=3)
        ax.text(peak_x, min(peak_y + y_limit * 0.05, y_limit * 0.96), f"{peak_x:.1f} Ma", ha="center", va="bottom", fontsize=6, color=base_color)
        ax.set_title(display_name_map.get(sample_name, sample_name).replace("_", " "), fontsize=7.5, pad=3, fontstyle="italic")
        ax.text(0.97, 0.9, f"n={curve_info['count']}", transform=ax.transAxes, ha="right", va="top", fontsize=6, color="#444444")
        ax.set_xlim(MIN_AGE, MAX_AGE)
        ax.set_ylim(0, y_limit)
        ax.grid(axis="y", color="#dddddd", linewidth=0.5, alpha=0.7)
        ax.tick_params(labelsize=6, length=2)
    for ax in axes[panel_count:]:
        ax.axis("off")
    fig.suptitle("Per-sample Duplication-Age Curves", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(output_dir / "combined_sample_line_panels.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_plot_combined(ages_dir: Path, clade_file: Path, tree_file: Path, output_dir: Path, mapping_file: Path | None = None) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_data = load_all_samples(ages_dir)
    if not sample_data:
        raise ValueError("当前目录没有可用于作图的 .pep 年龄文件，或数据点太少。")
    display_name_map = load_display_name_map(mapping_file)
    tree, leaf_order = read_tree_leaf_order(tree_file)
    sample_data, missing = reorder_sample_data_by_tree(sample_data, leaf_order)
    if missing:
        raise ValueError(f"树中存在没有年龄数据的样本: {missing[:5]}")
    clade_groups = read_clade_groups(clade_file)
    sample_to_clade = build_sample_to_clade(clade_groups)
    all_ages = np.concatenate([ages for _, ages in sample_data])
    bin_edges = np.linspace(MIN_AGE, MAX_AGE, BIN_COUNT + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    clade_densities = compute_clade_densities(sample_data, sample_to_clade, bin_edges)
    heatmap, sample_names, sample_sizes, peak_ages, ordered_clades = build_heatmap(sample_data, bin_edges, sample_to_clade)
    rgb_heatmap = build_colored_heatmap(heatmap, ordered_clades)

    fig = plt.figure(figsize=(20.5, max(8, len(sample_names) * 0.22 + 3)))
    gs = fig.add_gridspec(nrows=2, ncols=2, height_ratios=[1.2, 4], width_ratios=[4.8, 5.2], hspace=0.08, wspace=0.02)
    ax_top_left = fig.add_subplot(gs[0, 0]); ax_top_left.axis("off")
    ax_top = fig.add_subplot(gs[0, 1])
    ax_top.bar(bin_centers, np.histogram(all_ages, bins=bin_edges, density=True)[0], width=bin_edges[1]-bin_edges[0], color="#d9d9d9", edgecolor="#666666", linewidth=0.4, alpha=0.55)
    for clade_name, density in clade_densities.items():
        ax_top.plot(bin_centers, smooth_curve(density), color=CLADE_COLORS.get(clade_name, "#666666"), linewidth=2.2, label=clade_name)
    ax_top.set_xlim(MIN_AGE, MAX_AGE)
    ax_top.set_ylabel("Density")
    ax_top.set_title("Combined Duplication-Age Distribution")
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_top.legend(frameon=False, fontsize=9, ncol=3)

    ax_tree = fig.add_subplot(gs[1, 0])
    draw_tree_axis(ax_tree, tree, sample_names, display_name_map)
    ax_bottom = fig.add_subplot(gs[1, 1], sharex=ax_top, sharey=ax_tree)
    ax_bottom.imshow(rgb_heatmap, aspect="auto", interpolation="nearest", extent=[MIN_AGE, MAX_AGE, len(sample_names), 0])
    ax_bottom.set_xlabel("Age (Ma)")
    ax_bottom.set_xticks(bin_edges, minor=True)
    ax_bottom.set_yticks(np.arange(0, len(sample_names) + 1, 1), minor=True)
    ax_bottom.grid(which="minor", color="white", linewidth=0.25, alpha=0.7)
    ax_bottom.tick_params(which="minor", bottom=False, left=False)
    ax_bottom.tick_params(axis="y", left=False, labelleft=False)
    for spine in ax_bottom.spines.values():
        spine.set_visible(False)
    previous = None
    for index, clade_name in enumerate(ordered_clades):
        if previous is not None and clade_name != previous:
            ax_bottom.axhline(index, color="white", linewidth=1.2)
        previous = clade_name
    fig.savefig(output_dir / "combined_wgd_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    with open(output_dir / "combined_sample_stats.tsv", "w", encoding="utf-8") as handle:
        handle.write("sample\tclade\tcount\tpeak_age_ma\n")
        for sample_name, clade_name, sample_size, peak_age in zip(sample_names, ordered_clades, sample_sizes, peak_ages):
            handle.write(f"{sample_name}\t{clade_name}\t{sample_size}\t{peak_age:.4f}\n")

    plot_sample_line_panels(sample_data, sample_to_clade, display_name_map, bin_edges, output_dir)
