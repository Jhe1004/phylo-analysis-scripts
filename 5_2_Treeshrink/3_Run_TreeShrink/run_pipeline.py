import os
import glob
import re
import random
import subprocess
import shutil
from Bio import SeqIO

# ==============================================================================
# ---                           配置您的流程                           ---
# ==============================================================================

# --- 第1部分: 输入文件 ---
# 存放所有 RAxML 输出文件 (.tre) 和序列比对文件 (.fasta) 的文件夹
# "." 代表当前文件夹
INPUT_DIR = "."

# 序列比对文件的文件名模式 (使用通配符 *)
ALIGNMENT_PATTERN = "ortho*_cds_maffted.fas"

# RAxML 最佳树文件的文件名模式 (使用通配符 *)
# 请根据您实际的 RAxML 输出文件名进行修改
RAXML_TREE_PATTERN = "RAxML_bestTree.ortho*"


# --- 第2部分: 中间文件 ---
# 合并所有基因树后生成的临时文件名
COMBINED_TREE_FILE = "combined_trees.tre"


# --- 第3部分: TreeShrink 设置 ---
# 指向您电脑上的 run_treeshrink.py 脚本的路径
# 如果此流程脚本和 run_treeshrink.py 在同一个文件夹, 保持 "./run_treeshrink.py" 即可
TREESHRINK_SCRIPT_PATH = "./run_treeshrink.py"

# TreeShrink 运行时创建的输出文件夹名
TREESHRINK_OUTPUT_DIR = "result_treeshrink"

# TreeShrink 输出的、包含要删除物种列表的文件路径
# 通常是 "输出文件夹/output.txt"
TAXA_TO_REMOVE_FILE = os.path.join(TREESHRINK_OUTPUT_DIR, "output.txt")

# 为 TreeShrink 添加额外的命令行参数
# 例如, 使用 "per-species" 模式: TREESHRINK_EXTRA_ARGS = ["-a", "0.05"]
# 如果不需要额外参数, 保持为空列表: []
TREESHRINK_EXTRA_ARGS = []


# --- 第4部分: 最终输出 ---
# 存放过滤后序列文件的文件夹名
FINAL_ALIGNMENT_DIR = "shrunk_alignments"


# ==============================================================================
# ---                           脚本主逻辑 (通常无需修改)                      ---
# ==============================================================================

def extract_number_from_filename(filename):
    """从文件名 (如 'ortho76_cds_maffted.fasta') 中提取数字部分 (76)"""
    # 尝试多种正则表达式以提高兼容性
    patterns_to_try = [r'ortho(\d+)', r'(\d+).fasta', r'(\d+).tre']
    for pattern in patterns_to_try:
        match = re.search(pattern, os.path.basename(filename))
        if match:
            return int(match.group(1))
    return float('inf') # 如果找不到数字，则排在最后

def get_sorted_files(directory, pattern):
    """根据文件名中的数字查找并排序文件"""
    search_path = os.path.join(directory, pattern)
    files = glob.glob(search_path)
    if not files:
        print(f"警告: 在 '{directory}' 中未找到任何匹配 '{pattern}' 的文件。")
        return []
    return sorted(files, key=extract_number_from_filename)

def step1_combine_trees(input_dir, tree_pattern, output_file):
    """步骤一: 查找、排序并合并所有的基因树文件"""
    print("\n--- [步骤 1/3] 开始合并基因树 ---")
    tree_files = get_sorted_files(input_dir, tree_pattern)

    if not tree_files:
        print("错误: 未找到任何基因树文件，流程中止。")
        return False, None

    print(f"找到 {len(tree_files)} 个基因树文件，将按数字顺序合并到 '{output_file}'...")
    with open(output_file, "w") as outfile:
        for tree_file in tree_files:
            with open(tree_file, "r") as infile:
                content = infile.read().strip()
                if content: # 确保文件内容不为空
                    outfile.write(content + "\n")
    
    print("基因树合并完成！")
    return True, tree_files

def step2_run_treeshrink(script_path, combined_trees_file, output_dir, extra_args):
    """步骤二: 调用外部脚本 run_treeshrink.py"""
    print("\n--- [步骤 2/3] 开始运行 TreeShrink ---")
    if not os.path.exists(script_path):
        print(f"错误: TreeShrink 脚本在 '{script_path}' 未找到，流程中止。")
        return False
        
    # 构建命令行
    command = [
        "python", script_path,
        "-t", combined_trees_file,
        "-o", output_dir
    ] + extra_args

    print(f"将执行以下命令:\n{' '.join(command)}")
    
    try:
        # 执行命令并捕获输出
        result = subprocess.run(command, check=True, capture_output=True, text=True, encoding='utf-8')
        print("--- TreeShrink 输出 ---")
        print(result.stdout)
        print("-----------------------")
        print("TreeShrink 运行成功！")
        return True
    except subprocess.CalledProcessError as e:
        print("\n\033[91m错误: TreeShrink 运行失败！\033[0m")
        print("--- TreeShrink 错误信息 ---")
        print(e.stderr)
        print("---------------------------")
        return False

def step3_filter_alignments(alignment_dir, alignment_pattern, taxa_file, output_dir):
    """步骤三: 根据 TreeShrink 结果过滤序列文件"""
    print("\n--- [步骤 3/3] 开始根据 TreeShrink 结果过滤序列文件 ---")
    
    if not os.path.exists(taxa_file):
        print(f"错误: 找不到 TreeShrink 的输出文件 '{taxa_file}'，流程中止。")
        return

    # 这部分逻辑直接复用我们之前验证过的脚本
    alignment_files = get_sorted_files(alignment_dir, alignment_pattern)
    if not alignment_files:
        print("错误: 未找到任何序列比对文件，流程中止。")
        return

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    with open(taxa_file, 'r') as f:
        lines_to_process = f.readlines()

    if len(alignment_files) != len(lines_to_process):
        print(f"警告: 比对文件数量 ({len(alignment_files)}) 与 taxa 列表行数 ({len(lines_to_process)}) 不匹配!")

    verification_stats = {}
    print("开始处理每个比对文件...")
    for i, file_path in enumerate(alignment_files):
        if i >= len(lines_to_process):
            continue

        line = lines_to_process[i]
        taxa_to_remove = set(line.strip().split())
        
        original_records = list(SeqIO.parse(file_path, "fasta"))
        records_to_keep = [rec for rec in original_records if rec.id not in taxa_to_remove]
        
        output_filename = os.path.basename(file_path)
        output_path = os.path.join(output_dir, output_filename)
        SeqIO.write(records_to_keep, output_path, "fasta")

        verification_stats[file_path] = {
            'original_count': len(original_records),
            'removed_count_in_list': len(taxa_to_remove),
            'actually_removed': len(original_records) - len(records_to_keep),
            'actual_final_count': len(records_to_keep)
        }

    print("所有比对文件处理完毕！")
    verify_filter_results(verification_stats)


def verify_filter_results(stats, num_to_check=5):
    """随机抽查过滤结果并打印报告"""
    print("\n--- 验证过滤结果 ---")
    if not stats:
        print("没有可验证的文件。")
        return

    files_to_check = random.sample(list(stats.keys()), min(num_to_check, len(stats)))
    all_passed = True
    for file_path in files_to_check:
        stat = stats[file_path]
        filename = os.path.basename(file_path)
        
        print(f"\n--- 验证报告: {filename} ---")
        print(f"  原始序列数: {stat['original_count']}")
        print(f"  列表欲删除数: {stat['removed_count_in_list']}")
        print(f"  实际删除数: {stat['actually_removed']}")
        print(f"  最终序列数: {stat['actual_final_count']}")
        
        if stat['original_count'] - stat['actually_removed'] == stat['actual_final_count']:
            print("  状态: \033[92m成功\033[0m")
        else:
            print("  状态: \033[91m失败\033[0m")
            all_passed = False
            
    print("\n--- 验证总结 ---")
    if all_passed:
        print("\033[92m所有抽查文件均通过验证！\033[0m")
    else:
        print("\033[91m警告：部分抽查文件未能通过验证。\033[0m")

if __name__ == "__main__":
    print("="*50)
    print("      启动 RAxML -> TreeShrink -> Filter 自动化流程")
    print("="*50)

    # 步骤 1
    success, sorted_trees = step1_combine_trees(INPUT_DIR, RAXML_TREE_PATTERN, COMBINED_TREE_FILE)
    if not success:
        exit(1)

    # 步骤 2
    success = step2_run_treeshrink(TREESHRINK_SCRIPT_PATH, COMBINED_TREE_FILE, TREESHRINK_OUTPUT_DIR, TREESHRINK_EXTRA_ARGS)
    if not success:
        exit(1)

    # 步骤 3
    step3_filter_alignments(INPUT_DIR, ALIGNMENT_PATTERN, TAXA_TO_REMOVE_FILE, FINAL_ALIGNMENT_DIR)

    print("\n\033[92m====================================================\033[0m")
    print(f"\n\033[92m      流程全部执行完毕！\033[0m")
    print(f"\n\033[92m      最终过滤后的比对文件已保存至 '{FINAL_ALIGNMENT_DIR}' 文件夹。\033[0m")
    print("\n\033[92m====================================================\033[0m")
