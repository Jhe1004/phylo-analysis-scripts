#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import subprocess
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEPENDENCIES_DIR = SCRIPT_DIR / "dependencies"


def resolve_conda_bin(env_name: str) -> Path | None:
    for prefix in (
        Path.home() / "miniconda3" / "envs" / env_name / "bin",
        Path.home() / "mambaforge" / "envs" / env_name / "bin",
        Path.home() / "anaconda3" / "envs" / env_name / "bin",
        Path("/home/hejian/anaconda3/envs") / env_name / "bin",
    ):
        if prefix.exists():
            return prefix
    return None


def build_runtime_env(conda_env_name: str) -> dict[str, str]:
    env = os.environ.copy()
    path_parts = []

    local_bin = DEPENDENCIES_DIR / "bin"
    if local_bin.exists():
        path_parts.append(str(local_bin))

    conda_bin = resolve_conda_bin(conda_env_name)
    if conda_bin is not None:
        path_parts.append(str(conda_bin))

    path_parts.append(env.get("PATH", ""))
    env["PATH"] = ":".join(part for part in path_parts if part)
    return env


def ensure_exists(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"未找到{description}: {path}")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def run_shell_command(command: str, cwd: Path, env: dict[str, str]) -> None:
    print(f"\n[执行] {command}")
    subprocess.run(command, cwd=str(cwd), env=env, shell=True, check=True)
