## 使用步骤

1. 在Geneious软件中安装MAFFT插件，并将export_annotation.geneiousWorkfolw; extract_annotation.geneiousWorkfolw; extract_annotation.geneiousWorkfolw; 三个workfolw导入Geneious软件中.

2. 将所有叶绿体基因组的gb文件导入Geneious中。

3. 在Geneious中选中所有文件，Workflows中点击extract。该步骤会将叶绿体中所有的功能区提取出来，每个生成的文件中包含一个叶绿体的所有功能基因，并分行显示。

4. 选中所有生成的文件，Workflows中点击export。该步骤会将这些刚刚生成的文件导出到桌面的Batch Export文件夹中（文件夹会自动建立）。

5. 将Batch Export文件夹中的所有文件拷入本脚本所在的文件夹。

6. 运行脚本Formatting_Annotation.py 该脚本会将大多数的rRNA注释修改为普通的rRNA注释，并且会删掉所有的tRNA，因为tRNA注释信息不同的序列注释方法差异很大，难以统一，强行排序会引入很多的非同源错误。并且tRNA总长度只有两到三千，长度较短，删除后对最终结果影响很小。

7. 运行脚本Combine_Uniform_Gene.py 该脚本会将同名基因中的一个基因删掉，因为叶绿体中存在反向重复区，反向重复区基因序列完全一样，重复计算会加权该区域的比重。

8. 运行脚本Transform.py 该脚本会将所有文件中相同名称的基因序列提取出来，形成一个以基因名称命名的文件。需要注意，1. 脚本第一行写入了常规叶绿体中包含的基因名称，可以检查并且添加新的基因名称。 2. 该脚本还会生成一个response.txt文件，该文件中包含有各个基因文件中包含有哪些物种的基因情况。

9. 将这些结果导入Geneious软件中，选中所有文件，Workflows中点击Mafft。软件会自动将所有基因文件进行排序。

10. 排序的结果直接使用geneious软件的Tools - Concatenate Sequences or Alignments功能进行连接。

11. 使用RAxML软件进行建树。