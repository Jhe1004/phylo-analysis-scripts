#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import os

from ete3 import Tree


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"

TREE_FILE_PREFIX = "RAxML_bestTree."
OUTGROUP_FILE = "outgroupspecies.txt"
PROCESS_COUNT = max(1, min(8, os.cpu_count() or 1))
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_species_list(file_path: Path) -> set[str]:
    return {line.strip() for line in file_path.read_text(encoding="utf-8").splitlines() if line.strip()}


def get_non_outgroup_genes(gene_tree: Tree, outgroup_species: set[str]) -> list[str]:
    return [leaf for leaf in gene_tree.get_leaf_names() if leaf.split("++")[0] not in outgroup_species]


def get_outgroup_genes(gene_tree: Tree, outgroup_species: set[str]) -> list[str]:
    return [leaf for leaf in gene_tree.get_leaf_names() if leaf.split("++")[0] in outgroup_species]


def try_drop_one_intruding_outgroup(gene_tree: Tree, outgroup_species: set[str], non_outgroup_genes: list[str]):
    outgroup_genes = get_outgroup_genes(gene_tree, outgroup_species)
    if len(outgroup_genes) < 3:
        return None, None
    ancestor = gene_tree.get_common_ancestor(non_outgroup_genes)
    intruding = [leaf for leaf in ancestor.get_leaf_names() if leaf.split("++")[0] in outgroup_species]
    if len(intruding) != 1:
        return None, None
    pruned_tree = gene_tree.copy(method="deepcopy")
    pruned_tree.prune([leaf for leaf in pruned_tree.get_leaf_names() if leaf != intruding[0]], preserve_branch_length=True)
    pruned_non_outgroup = get_non_outgroup_genes(pruned_tree, outgroup_species)
    if not pruned_non_outgroup:
        return None, None
    pruned_ancestor = pruned_tree.get_common_ancestor(pruned_non_outgroup)
    if len(pruned_ancestor.get_leaf_names()) != len(pruned_non_outgroup):
        return None, None
    pruned_tree.set_outgroup(pruned_ancestor)
    return pruned_tree, intruding[0]


def reroot_tree(tree_path: Path, output_dir: Path, outgroup_species: set[str]) -> tuple[str, str, int, str]:
    try:
        gene_tree = Tree(str(tree_path))
    except Exception as exc:
        return tree_path.name, f"parse_error: {exc}", 0, "NA"

    non_outgroup_genes = get_non_outgroup_genes(gene_tree, outgroup_species)
    if not non_outgroup_genes:
        return tree_path.name, "no_non_outgroup_gene", 0, "NA"
    if len(non_outgroup_genes) == len(gene_tree.get_leaf_names()):
        return tree_path.name, "all_non_outgroup_no_outgroup", len(non_outgroup_genes), "NA"

    ancestor = gene_tree.get_common_ancestor(non_outgroup_genes)
    if len(ancestor.get_leaf_names()) != len(non_outgroup_genes):
        rescued_tree, deleted_leaf = try_drop_one_intruding_outgroup(gene_tree, outgroup_species, non_outgroup_genes)
        if rescued_tree is None:
            return tree_path.name, "non_outgroup_not_monophyletic", len(non_outgroup_genes), "NA"
        rescued_tree.write(format=1, outfile=str(output_dir / f"{tree_path.name}_r.tree"))
        return tree_path.name, "success_after_drop_one_outgroup", len(non_outgroup_genes), deleted_leaf

    gene_tree.set_outgroup(ancestor)
    gene_tree.write(format=1, outfile=str(output_dir / f"{tree_path.name}_r.tree"))
    return tree_path.name, "success", len(non_outgroup_genes), "NA"


def reroot_tree_task(task: tuple[Path, Path, set[str]]) -> tuple[str, str, int, str]:
    return reroot_tree(*task)


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    rooted_dir = output_dir / "rooted_gene_trees"
    ensure_directory(rooted_dir)

    outgroup_species = read_species_list(input_dir / OUTGROUP_FILE)
    tree_files = sorted(path for path in input_dir.iterdir() if path.is_file() and path.name.startswith(TREE_FILE_PREFIX))
    if not tree_files:
        raise FileNotFoundError(f"未找到前缀为 {TREE_FILE_PREFIX} 的树文件")

    tasks = [(tree_path, rooted_dir, outgroup_species) for tree_path in tree_files]
    if PROCESS_COUNT > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=PROCESS_COUNT) as executor:
            results = list(executor.map(reroot_tree_task, tasks))
    else:
        results = [reroot_tree_task(task) for task in tasks]

    counter = Counter(status for _, status, _, _ in results)
    log_file = output_dir / "reroot_results.tsv"
    with open(log_file, "w", encoding="utf-8") as handle:
        handle.write("tree\treason\tnon_outgroup_genes\tdeleted_leaf\n")
        for tree_name, reason, count, deleted_leaf in results:
            handle.write(f"{tree_name}\t{reason}\t{count}\t{deleted_leaf}\n")

    summary_file = output_dir / "reroot_summary.txt"
    with open(summary_file, "w", encoding="utf-8") as handle:
        for reason, count in sorted(counter.items()):
            handle.write(f"{reason}: {count}\n")

    print(f"完成树数: {len(results)}")
    print(f"成功输出目录: {rooted_dir}")


if __name__ == "__main__":
    main()
