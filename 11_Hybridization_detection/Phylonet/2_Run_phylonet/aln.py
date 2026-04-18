import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


# ==================== 配置区（直接修改这里）====================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/Phylonet/2_Run_phylonet/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/Phylonet/2_Run_phylonet/output"
CONDA_ENV_NAME = "trinity_env"

INPUT_EXTENSIONS = [".nex"]
PROCESS_COUNT = 8

PHYLONET_JAR = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/Phylonet/2_Run_phylonet/PhyloNet_3.8.0.jar"
JAVA_EXECUTABLE_NAME = "java"
# ============================================================


LOG_FILE = "run_phylonet.log"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY
DEPENDENCIES_DIR = SCRIPT_DIR / "dependencies"


def find_conda_env_root(env_name: str) -> Path | None:
    env_path = Path.home() / ".conda" / "envs" / env_name
    return env_path if env_path.exists() else None


def resolve_executable(executable_name: str, conda_env_name: str) -> str:
    local_candidate = DEPENDENCIES_DIR / "bin" / executable_name
    if local_candidate.exists():
        return str(local_candidate)

    env_root = find_conda_env_root(conda_env_name)
    if env_root:
        env_candidate = env_root / "bin" / executable_name
        if env_candidate.exists():
            return str(env_candidate)

    system_candidate = shutil.which(executable_name)
    if system_candidate:
        return system_candidate

    raise FileNotFoundError(f"未找到可执行程序 {executable_name}。")


def resolve_phylonet_jar(jar_name: str, conda_env_name: str) -> str:
    search_candidates = [
        DEPENDENCIES_DIR / "bin" / jar_name,
        DEPENDENCIES_DIR / jar_name,
        SCRIPT_DIR / jar_name,
    ]

    env_root = find_conda_env_root(conda_env_name)
    if env_root:
        search_candidates.extend(
            [
                env_root / "bin" / jar_name,
                env_root / "share" / jar_name,
                env_root / "lib" / jar_name,
                env_root / "libexec" / jar_name,
            ]
        )

    for candidate in search_candidates:
        if candidate.exists():
            return str(candidate)

    raise FileNotFoundError(
        "未找到 PhyloNet jar 文件。请将 jar 放到当前目录的 dependencies/bin/ 中，"
        f"或修改 PHYLONET_JAR 指向正确文件。当前目标文件名: {jar_name}"
    )


def run_single_nexus(java_path: str, jar_path: str, nexus_path: Path, output_path: Path) -> tuple[str, int]:
    with output_path.open("w", encoding="utf-8") as output_handle:
        process = subprocess.run(
            [java_path, "-jar", jar_path, str(nexus_path)],
            stdout=output_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return nexus_path.name, process.returncode


def main() -> None:
    if not INPUT_DIR.exists():
        raise FileNotFoundError(f"未找到输入目录: {INPUT_DIR}")

    nexus_files = sorted(
        path for path in INPUT_DIR.iterdir() if path.is_file() and path.suffix in INPUT_EXTENSIONS
    )
    if not nexus_files:
        raise FileNotFoundError(f"在 {INPUT_DIR} 中未找到 {INPUT_EXTENSIONS} 文件。")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    java_path = resolve_executable(JAVA_EXECUTABLE_NAME, CONDA_ENV_NAME)
    jar_path = resolve_phylonet_jar(PHYLONET_JAR, CONDA_ENV_NAME)
    log_path = OUTPUT_DIR / LOG_FILE

    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"JAVA: {java_path}\n")
        log_handle.write(f"PHYlONET_JAR: {jar_path}\n")
        log_handle.write(f"PROCESS_COUNT: {PROCESS_COUNT}\n")

        with ProcessPoolExecutor(max_workers=PROCESS_COUNT) as executor:
            future_map = {}
            for nexus_path in nexus_files:
                output_path = OUTPUT_DIR / f"{nexus_path.stem}.txt"
                future = executor.submit(run_single_nexus, java_path, jar_path, nexus_path, output_path)
                future_map[future] = nexus_path.name

            failures = []
            for future in as_completed(future_map):
                file_name, return_code = future.result()
                log_handle.write(f"{file_name}\t{return_code}\n")
                if return_code != 0:
                    failures.append(file_name)

    if failures:
        raise RuntimeError(f"以下输入文件运行失败: {', '.join(failures)}")

    print(f"PhyloNet 运行完成，共处理 {len(nexus_files)} 个输入文件。")
    print(f"日志文件: {log_path}")


if __name__ == "__main__":
    main()
