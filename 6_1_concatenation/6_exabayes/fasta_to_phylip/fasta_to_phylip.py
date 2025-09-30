import argparse
import re
import sys
from Bio import SeqIO

def main():
    """
    Main function to parse arguments and perform the conversion using BioPython.
    This version creates a stricter PHYLIP format for better compatibility.
    """
    parser = argparse.ArgumentParser(
        description="使用 BioPython 将 FASTA 格式的比对文件转换为兼容性更好的交错式 PHYLIP 格式。",
        epilog="使用示例: python fasta_to_phylip_strict.py your_alignment.fasta converted_alignment.phy"
    )
    parser.add_argument("input_fasta", help="输入的 FASTA 文件路径。")
    parser.add_argument("output_phylip", help="输出的 PHYLIP 文件路径。")
    args = parser.parse_args()

    # 1. 使用 BioPython 读取 FASTA 文件
    try:
        records = list(SeqIO.parse(args.input_fasta, "fasta"))
    except FileNotFoundError:
        print(f"错误: 输入文件 '{args.input_fasta}' 未找到。", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"解析 FASTA 文件时出错: {e}", file=sys.stderr)
        sys.exit(1)

    if not records:
        print("错误: 输入文件中未找到任何序列。", file=sys.stderr)
        sys.exit(1)

    # 2. 清理序列名并检查序列长度
    seq_len = len(records[0].seq)

    for record in records:
        if len(record.seq) != seq_len:
            print(f"错误: 序列长度不一致。序列 '{record.id}' 的长度为 {len(record.seq)}, 应为 {seq_len}。", file=sys.stderr)
            sys.exit(1)
        
        # 清理序列名
        sanitized_id = re.sub(r'[^a-zA-Z0-9_]', '_', record.id)
        
        # --- 关键改动 1: 将序列名截断或填充至正好10个字符 ---
        record.id = sanitized_id[:10].ljust(10)
        record.description = ""

    num_seqs = len(records)
    
    # 3. 写入 PHYLIP 文件
    block_size = 60
    sub_block_size = 10

    try:
        with open(args.output_phylip, 'w') as f:
            f.write(f" {num_seqs} {seq_len}\n") # PHYLIP标准通常在数字前有一个空格

            for i in range(0, seq_len, block_size):
                if i > 0:
                    f.write("\n")
                
                for record in records:
                    seq_str = str(record.seq)
                    current_block = seq_str[i : i + block_size]
                    
                    # --- 关键改动 2: 不再插入空格，以提高兼容性 ---
                    # formatted_block_parts = [current_block[j:j+sub_block_size] for j in range(0, len(current_block), sub_block_size)]
                    # formatted_block = ' '.join(formatted_block_parts)
                    
                    # 写入序列名（正好10个字符）+ 两个空格 + 序列数据
                    f.write(f"{record.id}  {current_block}\n")
                    
        print(f"成功将 '{args.input_fasta}' 转换为兼容性更好的 PHYLIP 文件 '{args.output_phylip}'。")

    except IOError as e:
        print(f"错误: 无法写入输出文件 '{args.output_phylip}': {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()