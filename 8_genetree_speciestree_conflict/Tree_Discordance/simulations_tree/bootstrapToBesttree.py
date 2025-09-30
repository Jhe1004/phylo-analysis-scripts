import subprocess

def run_raxml(t_file, z_file):
    # 完整命令
    command = [
        "raxmlHPC",  # 替换为RAxML的实际安装路径
        "-f", "b",
        "-m", "PROTGAMMAILG",
        "-n", "output_bootstrap.tre",
        "-t", t_file,
        "-z", z_file
    ]

    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("RAxML command executed successfully.")
        print("Output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error during RAxML execution:", e)
        print("Error Output:", e.stderr)

# 直接指定文件路径
t_file = "dna.tree"  # 物种树
z_file = "simulate.trees"  # Bootstrap树


# 执行函数
run_raxml(t_file, z_file)
