import glob
import os
import re

def process_gb_file(filename):
    """
    处理单个GenBank文件。
    """
    print(f"--- 处理中: {filename} ---")
    
    # 1. 从文件名生成 new_locus_name
    base_name = os.path.basename(filename)
    new_locus_name, _ = os.path.splitext(base_name)
    
    # 规则：替换所有非字母数字字符为下划线 '_'
    # [^a-zA-Z0-9] 匹配任何不是字母(a-z, A-Z)或数字(0-9)的字符
    new_locus_name = re.sub(r'[^a-zA-Z0-9]', '_', new_locus_name)
    
    # 清理可能出现的连续下划线 (例如 "a--b" 变为 "a_b")
    new_locus_name = re.sub(r'_+', '_', new_locus_name)
    
    # 移除开头或结尾可能多余的下划线
    new_locus_name = new_locus_name.strip('_')
    
    print(f"  - 新的LOCUS名称: {new_locus_name}")

    # 2. 读写文件
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"  - 错误: 无法读取文件 {filename}。原因: {e}")
        return

    if not lines:
        print(f"  - 警告: 文件 {filename} 为空。")
        return

    # 3. 修改 LOCUS 行 (通常是第一行)
    if lines[0].startswith('LOCUS'):
        original_line = lines[0]
        
        # 按空白符分割行
        parts = original_line.split()
        
        if len(parts) > 1:
            old_name = parts[1]
            
            # 替换第二个元素 (即旧的LOCUS名称)
            parts[1] = new_locus_name
            
            # 重新组合，用单个空格分隔
            # 这会丢失原始的列对齐，但是能确保新名称被完整替换
            lines[0] = " ".join(parts) + '\n'
            
            # 4. 写回文件
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    f.writelines(lines)
                print(f"  - 成功: {filename} 已更新, 旧名称 '{old_name}' 被替换。")
            except Exception as e:
                print(f"  - 错误: 无法写回文件 {filename}。原因: {e}")
        else:
            print(f"  - 警告: {filename} 的LOCUS行格式不正确。跳过。")
    else:
        print(f"  - 警告: {filename} 的第一行不是 'LOCUS' 行。跳过。")

def main():
    # 查找当前文件夹中所有的 .gb 文件
    gb_files = glob.glob('*.gb')
    
    if not gb_files:
        print("在当前目录中未找到任何 .gb 文件。")
        return

    print(f"总共找到 {len(gb_files)} 个 .gb 文件。")
    
    for gb_file in gb_files:
        process_gb_file(gb_file)
    
    print("\n--- 所有文件处理完毕 ---")

if __name__ == "__main__":
    main()