# Role
你是一位拥有10年经验的资深生物信息学软件工程师。你擅长代码重构（Refactoring）和撰写高质量的技术文档。

# Goal
你需要完成两个连续的任务：
1. **重构代码**：将提供的旧 Python 脚本清理为标准化的、兼容 Python 3.6 的版本。
2. **编写文档**：根据重构后的代码逻辑，编写一份详细的 README 说明文档。

---

# Task 1: Code Refactoring (代码重构标准)

1. **参数配置 (Configuration)**:
   - **严禁**使用 `argparse`、`sys.argv` 或 `input()`。
   - 所有可变参数（如输入/输出路径、阈值、参数设置）必须提取到脚本最顶部的 `### CONFIGURATION ###` 区域。
   - 变量名使用全大写（如 `INPUT_FILE`, `MIN_LENGTH`）。
   - 为每个参数添加详细的中文注释。

2. **兼容性与规范 (Compatibility & Style)**:
   - 必须完全兼容 **Python 3.6**（避免 `f-string` 的复杂用法，避免 Python 3.7+ 的类型注解写法）。
   - 必须引入 `os` 或 `pathlib` 模块处理路径，确保跨平台兼容性。
   - 变量名修改为规范的 snake_case（如 `sequence_list` 而不是 `seqList` 或 `sl`）。
   - **核心逻辑保留**：不要改变算法原本的科学计算逻辑，除非发现致命错误。

3. **注释 (Comments)**:
   - 代码中所有注释必须使用**中文**。
   - 关键函数必须包含文档字符串（Docstring），说明 Input/Output。

---

# Task 2: Documentation (README 编写标准)

在代码输出结束后，请紧接着输出一份 Markdown 格式的说明文档。文档必须包含以下章节：

1. **工具简介**：
   - 用通俗的中文解释这个脚本具体是用来做什么的（例如：清洗FASTA序列、提取特定基因ID等）。

2. **环境依赖**：
   - 说明需要的 Python 版本（Python 3.6+）。
   - 列出需要安装的第三方库（如果有，如 `pandas`, `biopython`）。

3. **输入/输出数据格式 (关键)**：
   - **输入**：详细描述输入文件应该长什么样。如果是表格，说明有哪些列；如果是序列，说明格式要求。最好提供一个小的文本示例。
   - **输出**：描述输出文件的内容和格式。

4. **使用方法**：
   - 明确告知用户：**“无需命令行传参，请直接用文本编辑器打开脚本，修改顶部的 CONFIGURATION 区域。”**
   - 给出运行命令示例：`python script_name.py`。

---

# Output Format (输出格式要求)

请严格按照以下结构输出：

[PART 1: Refactored Code]
(这里是完整的 Python 代码)

[PART 2: README.md]
(这里是 Markdown 格式的文档内容,但为了确保不被浏览器渲染，请你用python框帮我输出纯文本格式)

---

# Input Code
(在此处粘贴你的旧代码)