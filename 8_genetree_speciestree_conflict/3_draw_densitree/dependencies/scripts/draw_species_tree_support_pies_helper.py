#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import toyplot
import toyplot.svg
import toytree


SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def load_name_mapping(mapping_file: Path | None) -> dict[str, str]:
    if mapping_file is None or not mapping_file.exists():
        return {}
    mapping = {}
    with open(mapping_file, "r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            current_name = row["current_sample_name"].strip()
            standard_name = row["standard_latin_name"].strip()
            if current_name and standard_name:
                mapping[current_name] = standard_name
    return mapping


def rename_tree_tips(tree, name_mapping: dict[str, str]):
    if not name_mapping:
        return tree
    for node in tree.treenode.traverse():
        if node.is_leaf() and node.name in name_mapping:
            node.name = name_mapping[node.name]
    return tree


def load_tip_order(tip_order_file: Path | None, name_mapping: dict[str, str]) -> list[str] | None:
    if tip_order_file is None or not tip_order_file.exists():
        return None
    with open(tip_order_file, "r", encoding="utf-8") as handle:
        raw = [line.strip() for line in handle if line.strip()]
    return [name_mapping.get(name, name) for name in raw[::-1]]


def load_supports(support_csv: Path) -> dict[int, float]:
    support_by_idx = {}
    with open(support_csv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            support_by_idx[int(row["node_idx"])] = float(row["support_percent"])
    return support_by_idx


def pie_wedge_path(cx: float, cy: float, radius: float, start_deg: float, end_deg: float, steps: int = 32) -> str:
    if end_deg < start_deg:
        start_deg, end_deg = end_deg, start_deg
    points = [(cx, cy)]
    for step in range(steps + 1):
        fraction = step / steps
        degree = start_deg + (end_deg - start_deg) * fraction
        radian = math.radians(degree)
        points.append((cx + radius * math.cos(radian), cy + radius * math.sin(radian)))
    commands = [f"M {points[0][0]:.3f} {points[0][1]:.3f}"]
    for x_value, y_value in points[1:]:
        commands.append(f"L {x_value:.3f} {y_value:.3f}")
    commands.append("Z")
    return " ".join(commands)


def compute_adaptive_radii(
    pie_centers: list[tuple[float, float, float]],
    pie_scale: float,
    max_radius: float,
) -> list[tuple[float, float, float, float]]:
    coordinates = [(cx, cy) for cx, cy, _ in pie_centers]
    adaptive = []
    for index, (cx, cy, support) in enumerate(pie_centers):
        nearest = None
        for other_index, (ox, oy) in enumerate(coordinates):
            if index == other_index:
                continue
            distance = math.hypot(cx - ox, cy - oy)
            if nearest is None or distance < nearest:
                nearest = distance
        radius = max_radius if nearest is None else min(max_radius, nearest * pie_scale)
        adaptive.append((cx, cy, support, radius))
    return adaptive


def build_pie_group(pie_centers: list[tuple[float, float, float, float]]) -> ET.Element:
    group = ET.Element(
        f"{{{SVG_NS}}}g",
        {"id": "node-support-pies", "style": "pointer-events:none"},
    )
    for cx, cy, support_percent, radius in pie_centers:
        support = max(0.0, min(100.0, support_percent)) / 100.0
        group.append(
            ET.Element(
                f"{{{SVG_NS}}}circle",
                {
                    "cx": f"{cx:.3f}",
                    "cy": f"{cy:.3f}",
                    "r": f"{radius:.3f}",
                    "fill": "#d73027",
                    "fill-opacity": "0.88",
                    "stroke": "#ffffff",
                    "stroke-width": "0.6",
                },
            )
        )
        if support > 0.0:
            group.append(
                ET.Element(
                    f"{{{SVG_NS}}}path",
                    {
                        "d": pie_wedge_path(cx, cy, radius, -90.0, -90.0 + 360.0 * support),
                        "fill": "#4575b4",
                        "fill-opacity": "0.96",
                        "stroke": "none",
                    },
                )
            )
    return group


def add_legend(svg_root: ET.Element) -> None:
    legend = ET.Element(
        f"{{{SVG_NS}}}g",
        {"id": "legend", "transform": "translate(42,42)", "style": "pointer-events:none"},
    )
    legend.append(
        ET.Element(
            f"{{{SVG_NS}}}rect",
            {
                "x": "0",
                "y": "0",
                "width": "170",
                "height": "52",
                "rx": "8",
                "ry": "8",
                "fill": "#ffffff",
                "fill-opacity": "0.85",
                "stroke": "#444444",
                "stroke-width": "0.8",
            },
        )
    )

    for cy_value, fill, label in [("18", "#4575b4", "支持"), ("37", "#d73027", "不支持")]:
        legend.append(
            ET.Element(
                f"{{{SVG_NS}}}circle",
                {
                    "cx": "18",
                    "cy": cy_value,
                    "r": "7",
                    "fill": fill,
                    "stroke": "#ffffff",
                    "stroke-width": "0.6",
                },
            )
        )
        text = ET.Element(
            f"{{{SVG_NS}}}text",
            {
                "x": "32",
                "y": "22" if label == "支持" else "41",
                "fill": "#222222",
                "font-size": "12px",
                "font-family": "sans-serif",
            },
        )
        text.text = label
        legend.append(text)

    svg_root.append(legend)


def run_draw_species_tree_support_pies(
    species_tree_path: Path,
    gene_trees_path: Path,
    support_csv_path: Path,
    output_path: Path,
    tip_order_path: Path | None = None,
    name_mapping_path: Path | None = None,
    width: int = 1000,
    height: int = 1200,
    pie_max_radius: float = 11.0,
    pie_scale: float = 0.30,
) -> None:
    species_tree = toytree.tree(str(species_tree_path))
    gene_trees = toytree.mtree(str(gene_trees_path))
    name_mapping = load_name_mapping(name_mapping_path)
    species_tree = rename_tree_tips(species_tree, name_mapping)
    gene_trees.treelist = [rename_tree_tips(tree, name_mapping) for tree in gene_trees]

    tip_order = load_tip_order(tip_order_path, name_mapping)
    support_by_idx = load_supports(support_csv_path)

    canvas = toyplot.Canvas(width=width, height=height)
    axes = canvas.cartesian(xlabel="Relative time (root to tip)")
    gene_trees.draw_cloud_tree(
        axes=axes,
        fixed_order=tip_order,
        edge_style={
            "stroke": "#9a9a9a",
            "stroke-opacity": 0.16,
            "stroke-width": 1.0,
        },
    )
    _, _, mark = species_tree.draw(
        axes=axes,
        fixed_order=tip_order,
        edge_type="c",
        edge_style={
            "stroke": "black",
            "stroke-width": 1.4,
        },
        tip_labels_align=True,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    toyplot.svg.render(canvas, str(output_path))

    pie_centers = []
    for node in species_tree.treenode.traverse():
        if node.is_leaf() or node.idx not in support_by_idx:
            continue
        x_data, y_data = mark.ntable[node.idx]
        x_px = float(axes.project("x", [x_data])[0])
        y_px = float(axes.project("y", [y_data])[0])
        pie_centers.append((x_px, y_px, support_by_idx[node.idx]))

    adaptive_pies = compute_adaptive_radii(pie_centers, pie_scale=pie_scale, max_radius=pie_max_radius)
    tree = ET.parse(output_path)
    root = tree.getroot()
    root.append(build_pie_group(adaptive_pies))
    add_legend(root)
    tree.write(output_path, encoding="utf-8", xml_declaration=True)
    print(f"已写出支持率饼图: {output_path}")
