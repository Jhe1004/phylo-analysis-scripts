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
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/4_reroot_gene_family_trees/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/4_reroot_gene_family_trees/output"

SPECIES_TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/4_reroot_gene_family_trees/input/species_tree.newick"
INGROUP_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/4_reroot_gene_family_trees/input/ingroupspecies.txt"
GENE_TREE_SUBDIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/4_reroot_gene_family_trees/input/gene_trees"

GENE_TREE_PREFIX = "RAxML_bipartitions"
PROCESS_COUNT = max(1, min(8, os.cpu_count() or 1))
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_ingroup_species(file_path: Path) -> list[str]:
    return [line.strip() for line in file_path.read_text(encoding="utf-8").splitlines() if line.strip()]


def get_file_list(tree_dir: Path) -> list[Path]:
    return sorted(path for path in tree_dir.iterdir() if path.is_file() and path.name.startswith(GENE_TREE_PREFIX))


def get_outgroup_sort_list(species_tree_file: Path, ingroup_species: list[str]) -> list[list[str]]:
    outgroup_sort_list: list[list[str]] = []
    species_tree = Tree(str(species_tree_file))

    def walk(node) -> None:
        if node.is_leaf():
            return
        if len(ingroup_species) == len(node.get_leaf_names()):
            return
        if ingroup_species[0] in node.children[0].get_leaf_names():
            outgroup_sort_list.append(node.children[1].get_leaf_names())
            walk(node.children[0])
        else:
            outgroup_sort_list.append(node.children[0].get_leaf_names())
            walk(node.children[1])

    walk(species_tree)
    return outgroup_sort_list


def get_genes_in_species(hypothesis_root: list[str], gene_tree: Tree) -> list[str]:
    result = []
    for gene_name in gene_tree.get_leaf_names():
        species_name = gene_name.split("++")[0]
        if species_name in hypothesis_root:
            result.append(gene_name)
    return result


def reroot_single_tree(task: tuple[Path, Path, list[str], list[list[str]], Path]) -> tuple[str, str]:
    tree_path, species_tree_file, ingroup_species, outgroup_sort_list, output_dir = task
    try:
        _ = Tree(str(species_tree_file))
        gene_tree = Tree(str(tree_path))
    except Exception as exc:
        return tree_path.name, f"parse_error: {exc}"

    gene_tree_ingroup = [leaf for leaf in gene_tree.get_leaf_names() if leaf.split('++')[0] in ingroup_species]
    if not gene_tree_ingroup:
        return tree_path.name, "no_ingroup_gene"

    for hypothesis_root in outgroup_sort_list:
        candidate_genes = get_genes_in_species(hypothesis_root, gene_tree)
        if not candidate_genes:
            continue
        try:
            if len(candidate_genes) == 1:
                gene_tree.set_outgroup(candidate_genes[0])
            else:
                gene_tree.set_outgroup(gene_tree_ingroup[0])
                ancestor = gene_tree.get_common_ancestor(candidate_genes)
                if len(ancestor.get_leaf_names()) == len(candidate_genes):
                    gene_tree.set_outgroup(ancestor)
                else:
                    root_gene = candidate_genes[0]
                    gene_tree.set_outgroup(root_gene)
                    ancestor_len = len(gene_tree.get_common_ancestor(gene_tree_ingroup).get_leaf_names())
                    for candidate_gene in candidate_genes:
                        gene_tree.set_outgroup(candidate_gene)
                        current_len = len(gene_tree.get_common_ancestor(gene_tree_ingroup).get_leaf_names())
                        if current_len < ancestor_len:
                            ancestor_len = current_len
                            root_gene = candidate_gene
                    gene_tree.set_outgroup(root_gene)
            output_file = output_dir / f"{tree_path.name}_r.tree"
            gene_tree.write(format=1, outfile=str(output_file))
            return tree_path.name, "success"
        except Exception:
            continue
    return tree_path.name, "no_valid_root_found"


def reroot_single_tree_task(task: tuple[Path, Path, list[str], list[list[str]], Path]) -> tuple[str, str]:
    return reroot_single_tree(task)


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    rooted_dir = output_dir / "rooted_gene_trees"
    ensure_directory(rooted_dir)

    species_tree_file = input_dir / SPECIES_TREE_FILE
    ingroup_file = input_dir / INGROUP_FILE
    tree_dir = input_dir / GENE_TREE_SUBDIRECTORY
    if not species_tree_file.exists():
        raise FileNotFoundError(f"未找到物种树文件: {species_tree_file}")
    if not ingroup_file.exists():
        raise FileNotFoundError(f"未找到内类群文件: {ingroup_file}")

    ingroup_species = read_ingroup_species(ingroup_file)
    outgroup_sort_list = get_outgroup_sort_list(species_tree_file, ingroup_species)
    tree_files = get_file_list(tree_dir)
    if not tree_files:
        raise FileNotFoundError(f"未在 {tree_dir} 中找到前缀为 {GENE_TREE_PREFIX} 的基因树")

    tasks = [(tree_path, species_tree_file, ingroup_species, outgroup_sort_list, rooted_dir) for tree_path in tree_files]
    if PROCESS_COUNT > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=PROCESS_COUNT) as executor:
            results = list(executor.map(reroot_single_tree_task, tasks))
    else:
        results = [reroot_single_tree_task(task) for task in tasks]

    counter = Counter(status for _, status in results)
    log_file = output_dir / "reroot_results.tsv"
    with open(log_file, "w", encoding="utf-8") as handle:
        handle.write("tree\tstatus\n")
        for tree_name, status in results:
            handle.write(f"{tree_name}\t{status}\n")

    summary_file = output_dir / "reroot_summary.txt"
    with open(summary_file, "w", encoding="utf-8") as handle:
        for status, count in sorted(counter.items()):
            handle.write(f"{status}: {count}\n")

    print(f"完成树数: {len(results)}")
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
