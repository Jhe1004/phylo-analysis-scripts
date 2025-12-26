import os
import hashlib
import sys

def get_file_md5(filename, block_size=65536):
    """
    分块计算大文件的MD5值，避免一次性读入内存。
    """
    md5_hash = hashlib.md5()
    try:
        with open(filename, "rb") as f:
            # 循环读取文件块并更新hash
            for chunk in iter(lambda: f.read(block_size), b""):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()
    except IOError:
        # 文件不存在或无法读取
        return None
    except Exception as e:
        print(f"计算MD5时发生未知错误 ({filename}): {e}", file=sys.stderr)
        return "CALCULATION_ERROR"

def batch_md5_from_list(input_list_file, output_tsv_file):
    """
    从输入的文件列表（每行一个文件名）中读取文件名，
    按顺序计算MD5，并输出为TSV文件。
    """
    
    print(f"开始读取文件列表: {input_list_file}")
    
    # 1. 读取文件名列表
    try:
        with open(input_list_file, "r", encoding="utf-8") as f:
            # 读取所有行并去除每行末尾的换行符
            filenames = [line.strip() for line in f.readlines()]
            # 过滤掉空行
            filenames = [name for name in filenames if name]
            
    except FileNotFoundError:
        print(f"错误: 输入文件 '{input_list_file}' 未找到。", file=sys.stderr)
        print("请确保该文件与脚本在同一目录中。")
        return
    except Exception as e:
        print(f"读取文件列表时出错: {e}", file=sys.stderr)
        return

    if not filenames:
        print("文件列表为空，程序退出。")
        return

    print(f"从列表中读取了 {len(filenames)} 个文件名。开始计算MD5...")

    # 2. 遍历列表，计算MD5并写入TSV
    processed_count = 0
    try:
        with open(output_tsv_file, "w", encoding="utf-8") as f_out:
            # 严格按照 filenames 列表中的顺序
            for filename in filenames:
                print(f"  > 正在处理: {filename}")
                
                # 检查文件是否存在
                if not os.path.isfile(filename):
                    print(f"  ! 警告: 文件未找到: {filename}")
                    md5_sum = "FILE_NOT_FOUND"
                else:
                    md5_sum = get_file_md5(filename)
                    if md5_sum is None:
                        # 这种情况通常是 "Permission Denied"
                        md5_sum = "READ_ERROR"

                # 写入TSV（使用制表符 \t 分隔）
                f_out.write(f"{filename}\t{md5_sum}\n")
                processed_count += 1

        print(f"\n处理完成！")
        print(f"成功处理了 {processed_count} / {len(filenames)} 个文件。")
        print(f"结果已保存到: {output_tsv_file}")

    except IOError as e:
        print(f"错误: 无法写入输出文件 {output_tsv_file}: {e}", file=sys.stderr)

# --- 脚本主入口 ---
if __name__ == "__main__":
    # --- 配置 ---
    # 包含文件名的输入文件 (来自您上传的 )
    INPUT_FILE_LIST = "md5.txt" 
    # 您希望生成的TSV输出文件
    OUTPUT_TSV_FILE = "md5_results.tsv"
    # --- 结束配置 ---
    
    batch_md5_from_list(INPUT_FILE_LIST, OUTPUT_TSV_FILE)