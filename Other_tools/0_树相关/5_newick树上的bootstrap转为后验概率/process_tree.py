import re
import os

def convert_support_to_probability(input_file, output_file):
    """
    读取一个Newick格式的树文件，将分支支持率从0-100的范围
    转换到0-1的概率范围，并将结果写入一个新文件。

    Args:
        input_file (str): 输入的Newick文件名。
        output_file (str): 输出的Newick文件名。
    """
    try:
        with open(input_file, 'r') as f:
            newick_string = f.read()

        # 这个函数会被正则表达式的sub方法调用
        # 它会提取匹配到的支持率，除以100，然后返回替换后的字符串
        def replace_support(match):
            support_value = float(match.group(1))
            probability = support_value / 100.0
            # 使用f-string来格式化输出，确保小数的精度
            return f"){probability}"

        # 这个正则表达式会查找这样的模式：一个右括号")"，后面跟着数字（支持率），
        # 然后跟着一个冒号":"。这可以准确地定位到内部节点的支持率。
        # re.sub会找到所有匹配项，并用replace_support函数的返回值进行替换
        modified_newick_string = re.sub(r"\)([\d\.]+)", replace_support, newick_string)

        with open(output_file, 'w') as f:
            f.write(modified_newick_string)

        print(f"处理完成！")
        print(f"输入文件: {input_file}")
        print(f"输出文件: {output_file}")
        print("支持率已成功转换为后验概率。")

    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file}")
    except Exception as e:
        print(f"处理过程中发生错误: {e}")

# --- 主程序 ---
# 定义输入和输出文件名
input_filename = 'ExaBayes_topologies.run-0.its consensus.newick'
output_filename = 'ExaBayes_topologies.run-0.its consensus_prob.newick'

# 检查输入文件是否存在
if os.path.exists(input_filename):
    # 运行转换函数
    convert_support_to_probability(input_filename, output_filename)
else:
    print(f"错误: 输入文件 '{input_filename}' 不在当前目录中。")
    print("请确保脚本和newick文件在同一个文件夹下。")