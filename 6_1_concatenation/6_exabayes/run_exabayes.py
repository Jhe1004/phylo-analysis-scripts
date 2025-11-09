#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
一个用于启动 ExaBayes mpirun 作业的 Python 脚本。
它会自动生成一个随机种子来替代 $RANDOM。
"""

import subprocess
import random
import sys

# 1. 生成随机种子
# $RANDOM 在 bash 中通常是 0 到 32767 之间的一个整数
random_seed = random.randint(0, 32767)

# 2. 定义您的命令
# mpirun -np 240 exabayes -f complete.phy -q partitions.txt -n combine -s $RANDOM -c config.nex
# 我们将其拆分为一个列表，这是 subprocess 推荐的做法，更安全
command_list = [
    "mpirun",
    "-np", "40",
    "exabayes",
    "-f", "complete.phy",
    "-q", "partitions.txt",
    "-n", "combine",
    "-s", str(random_seed),  # 将随机种子转为字符串
    "-c", "config.nex"
]

# 3. 运行命令
print("--- 准备执行 ExaBayes 命令 ---")
print(f"Running: {' '.join(command_list)}")
print("---------------------------------\n")

try:
    # check=True 使得如果 mpirun 执行失败 (返回非0退出码)，Python 会抛出异常
    # 实时显示输出，因为我们没有捕获 stdout/stderr
    result = subprocess.run(command_list, check=True, text=True)
    
    print("\n---------------------------------")
    print("--- 命令执行成功 ---")

except FileNotFoundError:
    print(f"\n[错误] 命令 '{command_list[0]}' 未找到。")
    print("请确保 mpirun (OpenMPI/MPICH) 和 exabayes 已经正确安装并位于您的系统 PATH 路径中。")
except subprocess.CalledProcessError as e:
    # 如果 mpirun 返回了错误
    print(f"\n[错误] 命令执行失败，退出码: {e.returncode}")
    print("请检查 ExaBayes 的输出日志或错误信息。")
except KeyboardInterrupt:
    # 捕获用户按 Ctrl+C 的操作
    print("\n\n[用户中断] 操作已取消。")
    sys.exit(1)
except Exception as e:
    # 捕获其他未知错误
    print(f"\n[未知错误] 发生异常: {e}")
    sys.exit(1)