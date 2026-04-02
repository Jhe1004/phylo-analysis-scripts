from pathlib import Path


# ==================== 配置区（直接修改这里）====================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"

INPUT_EXTENSIONS = [".txt"]
TREE_OUTPUT_SUBDIRECTORY = "tree"
PROBABILITY_OUTPUT_SUBDIRECTORY = "probability"
PROBABILITY_SUMMARY_FILE = "result.txt"
# ============================================================


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def extract_tree_line(result_file: Path) -> str | None:
    for line in result_file.read_text(encoding="utf-8", errors="ignore").splitlines():
        if "Visualize in Dendroscope" in line:
            return line.split("Visualize in Dendroscope", 1)[1].strip()
    return None


def extract_probability_line(result_file: Path) -> str | None:
    for line in result_file.read_text(encoding="utf-8", errors="ignore").splitlines():
        if "Total log probability" in line:
            return line.split("Total log probability", 1)[1].strip(" :\t")
    return None


def main() -> None:
    if not INPUT_DIR.exists():
        raise FileNotFoundError(f"未找到输入目录: {INPUT_DIR}")

    result_files = sorted(
        path for path in INPUT_DIR.iterdir() if path.is_file() and path.suffix in INPUT_EXTENSIONS
    )
    if not result_files:
        raise FileNotFoundError(f"在 {INPUT_DIR} 中未找到 {INPUT_EXTENSIONS} 结果文件。")

    tree_dir = OUTPUT_DIR / TREE_OUTPUT_SUBDIRECTORY
    probability_dir = OUTPUT_DIR / PROBABILITY_OUTPUT_SUBDIRECTORY
    tree_dir.mkdir(parents=True, exist_ok=True)
    probability_dir.mkdir(parents=True, exist_ok=True)

    summary_lines = []
    for result_file in result_files:
        tree_line = extract_tree_line(result_file)
        probability_line = extract_probability_line(result_file)

        if tree_line:
            (tree_dir / result_file.name).write_text(tree_line + "\n", encoding="utf-8")
        if probability_line:
            summary_lines.append(f"{result_file.name},{probability_line}")

    summary_path = probability_dir / PROBABILITY_SUMMARY_FILE
    summary_path.write_text(
        "\n".join(summary_lines) + ("\n" if summary_lines else ""),
        encoding="utf-8",
    )

    print(f"解析完成，共处理 {len(result_files)} 个结果文件。")
    print(f"树文件输出目录: {tree_dir}")
    print(f"概率汇总文件: {summary_path}")


if __name__ == "__main__":
    main()
