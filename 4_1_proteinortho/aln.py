import os
import glob

# 获取当前目录下所有.pep文件
pep_files = glob.glob("*.pep")

# 对每个pep文件进行清理
for pep_file in pep_files:
    # 创建临时文件
    temp_file = pep_file + '.tmp'
    
    # 读取原文件并处理
    with open(pep_file, 'r') as f_in, open(temp_file, 'w') as f_out:
        for line in f_in:
            # 如果是序列行（不是以>开头），移除*字符
            if not line.startswith('>'):
                line = line.replace('*', '')
            f_out.write(line)
    
    # 用处理后的文件替换原文件
    os.replace(temp_file, pep_file)
    print(f"已清理文件: {pep_file}")

# 运行proteinortho
os.system("proteinortho *.pep")