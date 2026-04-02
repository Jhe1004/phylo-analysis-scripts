import os
import glob
import logging
from datetime import datetime
from collections import Counter


def setup_logger(script_name):
    """配置日志记录器，同时输出到文件和控制台"""
    log_file = f"{script_name}.log"
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.DEBUG)
    
    logger.handlers.clear()
    
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    logger.info(f"日志文件: {log_file}")
    return logger


logger = setup_logger("5_check_non_atcg")

def is_fasta(filename):
    return filename.lower().endswith(('.fasta', '.fa', '.fna', '.fas'))

def count_non_atcg(filepath):
    """
    统计单个FASTA文件中的非ATCG字符。
    """
    total_bases = 0
    base_counts = Counter()
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>') or not line:
                    continue
                
                # 转换为大写以统一统计
                seq_upper = line.upper()
                total_bases += len(seq_upper)
                base_counts.update(seq_upper)
                
    except Exception as e:
        logger.error(f"读取文件 {os.path.basename(filepath)} 时出错: {e}")
        return None

    # 计算 ATCG 的总数
    atcg_count = base_counts['A'] + base_counts['T'] + base_counts['C'] + base_counts['G']
    # 非 ATCG 的总数 = 总碱基数 - ATCG数
    non_atcg_count = total_bases - atcg_count
    
    return total_bases, non_atcg_count

def main():
    # 获取当前目录下所有的fasta文件
    current_dir = os.getcwd()
    fasta_files = [f for f in os.listdir(current_dir) if os.path.isfile(f) and is_fasta(f)]
    fasta_files.sort()

    if not fasta_files:
        logger.info("当前文件夹中没有找到 .fasta / .fa / .fna / .fas 文件。")
        return

    table_data = []
    grand_total_bases = 0
    grand_total_non_atcg = 0

    logger.info(f"正在分析当前文件夹: {current_dir}\n")

    for file in fasta_files:
        result = count_non_atcg(file)
        if result:
            total, non_atcg = result
            grand_total_bases += total
            grand_total_non_atcg += non_atcg
            
            ratio = (non_atcg / total * 100) if total > 0 else 0
            table_data.append([file, total, non_atcg, ratio])

    # --- 打印表格 ---
    header = ["Filename", "Total Bases", "Non-ATCG", "Ratio (%)"]
    col_widths = [len(h) for h in header]
    for row in table_data:
        col_widths[0] = max(col_widths[0], len(str(row[0])))
        col_widths[1] = max(col_widths[1], len(f"{row[1]:,}"))
        col_widths[2] = max(col_widths[2], len(f"{row[2]:,}"))
        col_widths[3] = max(col_widths[3], len(f"{row[3]:.4f}"))

    row_format = " | ".join([f"{{:<{col_widths[0]}}}", f"{{:>{col_widths[1]}}}", f"{{:>{col_widths[2]}}}", f"{{:>{col_widths[3]}}}"])

    separator = "-" * (sum(col_widths) + 3 * 3)
    logger.info(separator)
    logger.info(row_format.format(*header))
    logger.info(separator)

    for row in table_data:
        logger.info(row_format.format(
            row[0], 
            f"{row[1]:,}", 
            f"{row[2]:,}", 
            f"{row[3]:.4f}"
        ))

    logger.info(separator)
    grand_ratio = (grand_total_non_atcg / grand_total_bases * 100) if grand_total_bases > 0 else 0
    logger.info(row_format.format(
        "TOTAL SUMMARY", 
        f"{grand_total_bases:,}", 
        f"{grand_total_non_atcg:,}", 
        f"{grand_ratio:.4f}"
    ))
    logger.info(separator)

if __name__ == '__main__':
    main()