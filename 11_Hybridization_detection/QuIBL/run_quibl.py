import shutil
import subprocess
import sys
from configparser import ConfigParser
from pathlib import Path


# ==================== 配置区（直接修改这里）====================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/QuIBL/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/QuIBL/output"
CONDA_ENV_NAME = "base"

TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/QuIBL/input/gene_trees.trees"
TOTAL_OUTGROUP = "Symphalangus.fasta.transdecoder.pep"
NUM_DISTRIBUTIONS = 2
LIKELIHOOD_THRESHOLD = 0.01
NUM_STEPS = 50
GRAD_ASCENT_SCALAR = 0.5
MULTIPROCESS = True
MAX_CORES = 20

PYTHON_EXECUTABLE_NAME = "python"
# ============================================================


OUTPUT_CSV_FILE = "quibl_results.csv"
RUNTIME_CONFIG_FILE = "quibl_runtime.ini"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY
DEPENDENCIES_DIR = SCRIPT_DIR / "dependencies"
UPSTREAM_SCRIPT = DEPENDENCIES_DIR / "upstream" / "QuIBL.py"


def find_conda_env_root(env_name: str) -> Path | None:
    env_path = Path.home() / ".conda" / "envs" / env_name
    return env_path if env_path.exists() else None


def resolve_python_executable(executable_name: str, conda_env_name: str) -> str:
    local_candidate = DEPENDENCIES_DIR / "bin" / executable_name
    if local_candidate.exists():
        return str(local_candidate)

    env_root = find_conda_env_root(conda_env_name)
    if env_root:
        env_candidate = env_root / "bin" / executable_name
        if env_candidate.exists():
            return str(env_candidate)

    current_python = Path(sys.executable)
    if current_python.exists():
        return str(current_python)

    system_candidate = shutil.which(executable_name)
    if system_candidate:
        return system_candidate

    python3_candidate = shutil.which("python3")
    if python3_candidate:
        return python3_candidate

    raise FileNotFoundError("未找到可用于运行 QuIBL 的 Python 解释器。")


def write_runtime_config(config_path: Path, tree_path: Path, output_csv_path: Path) -> None:
    config = ConfigParser()
    config["Input"] = {
        "treefile": str(tree_path),
        "numdistributions": str(NUM_DISTRIBUTIONS),
        "likelihoodthresh": str(LIKELIHOOD_THRESHOLD),
        "numsteps": str(NUM_STEPS),
        "gradascentscalar": str(GRAD_ASCENT_SCALAR),
        "totaloutgroup": TOTAL_OUTGROUP,
        "multiproc": str(MULTIPROCESS),
        "maxcores": str(MAX_CORES),
    }
    config["Output"] = {"OutputPath": str(output_csv_path)}
    with config_path.open("w", encoding="utf-8") as handle:
        config.write(handle)


def main() -> None:
    if not UPSTREAM_SCRIPT.exists():
        raise FileNotFoundError(f"未找到 QuIBL 上游脚本: {UPSTREAM_SCRIPT}")

    tree_path = INPUT_DIR / TREE_FILE
    if not tree_path.exists():
        raise FileNotFoundError(f"未找到输入树文件: {tree_path}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    runtime_config_path = OUTPUT_DIR / RUNTIME_CONFIG_FILE
    output_csv_path = OUTPUT_DIR / OUTPUT_CSV_FILE
    log_path = OUTPUT_DIR / "run_quibl.log"

    write_runtime_config(runtime_config_path, tree_path, output_csv_path)
    python_path = resolve_python_executable(PYTHON_EXECUTABLE_NAME, CONDA_ENV_NAME)

    command = [python_path, str(UPSTREAM_SCRIPT), str(runtime_config_path)]
    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"PYTHON: {python_path}\n")
        log_handle.write(f"UPSTREAM_SCRIPT: {UPSTREAM_SCRIPT}\n")
        log_handle.write(f"RUNTIME_CONFIG: {runtime_config_path}\n")
        process = subprocess.run(
            command,
            cwd=SCRIPT_DIR,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )

    if process.returncode != 0:
        raise RuntimeError(f"QuIBL 运行失败，请检查日志文件: {log_path}")

    print(f"QuIBL 运行完成，结果文件: {output_csv_path}")
    print(f"日志文件: {log_path}")


if __name__ == "__main__":
    main()
