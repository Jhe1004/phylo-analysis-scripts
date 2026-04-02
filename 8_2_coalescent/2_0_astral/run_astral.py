import os
import shutil
import subprocess
import sys


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
CONDA_ENV_NAME = "trinity_env"

INPUT_EXTENSION = ".trees"
ASTRAL_EXECUTABLE_NAME = "astral"
ASTRAL_JAR = ""
OUTGROUP = ""
ASTRAL_EXTRA_ARGS = ["-t", "3"]
ASTRAL_CPU_THREADS = int(os.environ.get("ASTRAL_CPU_THREADS", "80"))
DRY_RUN = False
# ============================================================


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def find_executable_in_conda_env(executable_name, env_name):
    try:
        result = subprocess.run(
            ["conda", "run", "-n", env_name, "which", executable_name],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=True,
        )
    except Exception:
        return None
    lines = result.stdout.strip().splitlines()
    if not lines:
        return None
    candidate = lines[-1].strip()
    if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    return None


def resolve_java():
    local_candidate = os.path.join(SCRIPT_DIR, "dependencies", "bin", "java")
    if os.path.isfile(local_candidate) and os.access(local_candidate, os.X_OK):
        return local_candidate
    conda_candidate = find_executable_in_conda_env("java", CONDA_ENV_NAME)
    if conda_candidate:
        return conda_candidate
    path_candidate = shutil.which("java")
    if path_candidate:
        return path_candidate
    raise FileNotFoundError("未找到可执行程序: java")


def resolve_astral_command():
    local_candidate = os.path.join(SCRIPT_DIR, "dependencies", "bin", ASTRAL_EXECUTABLE_NAME)
    if os.path.isfile(local_candidate) and os.access(local_candidate, os.X_OK):
        return [local_candidate]
    conda_candidate = find_executable_in_conda_env(ASTRAL_EXECUTABLE_NAME, CONDA_ENV_NAME)
    if conda_candidate:
        return [conda_candidate]
    path_candidate = shutil.which(ASTRAL_EXECUTABLE_NAME)
    if path_candidate:
        return [path_candidate]
    if ASTRAL_JAR:
        java_executable = resolve_java()
        return [java_executable, "-jar", ASTRAL_JAR]
    raise FileNotFoundError("没有找到可用的 ASTRAL。请配置 ASTRAL_EXECUTABLE_NAME 或 ASTRAL_JAR。")


def list_input_files(input_dir):
    return sorted(file_name for file_name in os.listdir(input_dir) if file_name.endswith(INPUT_EXTENSION))


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    os.makedirs(output_dir, exist_ok=True)
    if not os.path.isdir(input_dir):
        print(f"错误: 输入目录不存在: {input_dir}")
        sys.exit(1)
    input_files = list_input_files(input_dir)
    if not input_files:
        print(f"未找到后缀为 {INPUT_EXTENSION} 的输入文件。")
        sys.exit(0)
    try:
        command_prefix = resolve_astral_command()
    except FileNotFoundError as exc:
        print(f"错误: {exc}")
        sys.exit(1)
    for input_file_name in input_files:
        input_path = os.path.join(input_dir, input_file_name)
        output_file = os.path.join(output_dir, f"{os.path.splitext(input_file_name)[0]}.tree")
        command = command_prefix + ASTRAL_EXTRA_ARGS
        if ASTRAL_CPU_THREADS >= 2:
            command += ["-T", str(ASTRAL_CPU_THREADS), "-C"]
        if OUTGROUP:
            command += ["--outgroup", OUTGROUP]
        command += ["-i", input_path, "-o", output_file]
        print(f"[INFO] 运行: {input_file_name} -> {os.path.basename(output_file)}")
        if DRY_RUN:
            continue
        subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
