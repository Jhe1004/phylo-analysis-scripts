from Bio import Entrez
import time

# 设置你的邮箱（NCBI要求必须提供）
Entrez.email = "j.he930724@gmail.com"  # 替换为你的真实邮箱

# 输入你的序列号列表
accession_list = [
"OR801208","OR786464","OR801209","OR801211","OR801210","ON520699","OR801212","OR801213","OR801214","OR801215","OR801216","OR801219","OR801217","OR801218","OR801220","ON520701","ON411429","ON411431","OR801221","OR801222","OR801223","OR801224","PP998020","PP998021","PP998022","PP998023","PP998024","ON520702","PP998025","PP998026","PP998027","PP998028","PP998029","PP998030","PP998031","PP998032","ON520703","ON411438","PP998033","PP998034","ON520704","PP998035","PP998036","PP998037","PP998038","PP998039","PP998040","PP998041","ON520705","PP998042","PP998043","PP998044","PP998045","PP998046","PP998047","ON520700"
]

def download_genbank(accessions, file_format="gb"):
    """
    下载GenBank格式文件
    参数:
        accessions: 序列号列表
        file_format: 下载格式（默认gb即GenBank格式）
    """
    for acc in accessions:
        try:
            # 下载数据
            handle = Entrez.efetch(
                db="nucleotide",
                id=acc,
                rettype=file_format,
                retmode="text"
            )
            
            # 读取数据并保存
            with open(f"{acc}.{file_format}", "w") as f:
                f.write(handle.read())
            
            print(f"成功下载 {acc}.{file_format}")
            
            # 遵守NCBI的请求频率要求（每秒不超过3次）
            time.sleep(0.34)  # 设置适当延迟
            
        except Exception as e:
            print(f"下载 {acc} 失败: {str(e)}")

if __name__ == "__main__":
    download_genbank(accession_list)