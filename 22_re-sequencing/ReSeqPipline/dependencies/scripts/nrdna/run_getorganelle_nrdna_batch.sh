#!/usr/bin/env bash

set -euo pipefail

# 批量运行 GetOrganelle 组装 nrDNA（双端测序）
# 默认扫描当前目录下的 *.fq.gz 文件，并按 *_1 / *_2 自动配对
#
# 默认参数参考 GetOrganelle 官方常用 nrDNA 用法：
#   -F embplant_nr -R 10 -k 35,85,115
#
# 用法示例：
#   bash run_getorganelle_nrdna_batch.sh
#   bash run_getorganelle_nrdna_batch.sh 32
#   THREADS=40 OUTDIR=results bash run_getorganelle_nrdna_batch.sh

THREADS="${1:-${THREADS:-16}}"
OUTDIR="${OUTDIR:-getorganelle_nrdna_results}"
LOGDIR="${LOGDIR:-${OUTDIR}/logs}"
ENV_NAME="${ENV_NAME:-getorganelle}"
GETORG_CMD="${GETORG_CMD:-}"
GETORG_ENTRY="${GETORG_ENTRY:-}"
ORGANELLE_TYPE="${ORGANELLE_TYPE:-embplant_nr}"
KMER_LIST="${KMER_LIST:-35,85,115}"
MAX_ROUNDS="${MAX_ROUNDS:-10}"
EXTRA_ARGS="${EXTRA_ARGS:-}"

if ! command -v conda >/dev/null 2>&1; then
    echo "错误：未找到 conda，请先确保 conda 可用。"
    exit 1
fi

mkdir -p "${OUTDIR}" "${LOGDIR}"

resolve_getorg_cmd() {
    if [[ -n "${GETORG_CMD}" ]]; then
        printf '%s\n' "${GETORG_CMD}"
        return 0
    fi

    if conda run -n "${ENV_NAME}" get_organelle_from_reads.py -h >/dev/null 2>&1; then
        printf 'conda run -n %q' "${ENV_NAME}"
        return 0
    fi

    if conda run -n "${ENV_NAME}" get_organelle_from_reads -h >/dev/null 2>&1; then
        printf 'conda run -n %q' "${ENV_NAME}"
        return 0
    fi

    return 1
}

resolve_getorg_entry() {
    if [[ -n "${GETORG_ENTRY}" ]]; then
        printf '%s\n' "${GETORG_ENTRY}"
        return 0
    fi

    if conda run -n "${ENV_NAME}" get_organelle_from_reads.py -h >/dev/null 2>&1; then
        printf 'get_organelle_from_reads.py'
        return 0
    fi

    if conda run -n "${ENV_NAME}" get_organelle_from_reads -h >/dev/null 2>&1; then
        printf 'get_organelle_from_reads'
        return 0
    fi

    return 1
}

GETORG_CMD="$(resolve_getorg_cmd || true)"
GETORG_ENTRY="$(resolve_getorg_entry || true)"
if [[ -z "${GETORG_CMD}" || -z "${GETORG_ENTRY}" ]]; then
    echo "错误：在 conda 环境 ${ENV_NAME} 中没有找到 GetOrganelle 主程序。"
    echo "可手动指定，例如："
    echo "  GETORG_CMD='conda run -n ${ENV_NAME}' GETORG_ENTRY='get_organelle_from_reads.py' bash $0"
    exit 1
fi

echo "输出目录: ${OUTDIR}"
echo "日志目录: ${LOGDIR}"
echo "线程数: ${THREADS}"
echo "GetOrganelle 命令: ${GETORG_CMD} ${GETORG_ENTRY}"
echo "组装类型: ${ORGANELLE_TYPE}"
echo "k-mer: ${KMER_LIST}"
echo "最大迭代轮数: ${MAX_ROUNDS}"
echo

found=0
done_n=0
skip_n=0
fail_n=0

while IFS= read -r -d '' fq1; do
    found=1
    base_name="$(basename "${fq1}")"
    prefix="${base_name%_1.*fq.gz}"
    fq2="./${prefix}_2${base_name#${prefix}_1}"

    if [[ ! -f "${fq2}" ]]; then
        echo "未找到配对文件，跳过: ${fq1}"
        ((skip_n++))
        continue
    fi

    sample="${prefix}"
    sample_out="${OUTDIR}/${sample}"
    sample_log="${LOGDIR}/${sample}.log"

    if [[ -s "${sample_out}/extended_path_sequence.fasta" || -s "${sample_out}/final_assembly_graph.fastg" ]]; then
        echo "已存在结果，跳过: ${sample}"
        ((skip_n++))
        continue
    fi

    mkdir -p "${sample_out}"

    echo "============================================================"
    echo "开始样本: ${sample}"
    echo "R1: ${fq1}"
    echo "R2: ${fq2}"
    echo "输出: ${sample_out}"
    echo "日志: ${sample_log}"

    set +e
    bash -lc "${GETORG_CMD} ${GETORG_ENTRY@Q} -1 ${fq1@Q} -2 ${fq2@Q} -o ${sample_out@Q} -F ${ORGANELLE_TYPE@Q} -R ${MAX_ROUNDS@Q} -k ${KMER_LIST@Q} -t ${THREADS@Q} ${EXTRA_ARGS}" \
        >"${sample_log}" 2>&1
    status=$?
    set -e

    if [[ ${status} -eq 0 ]]; then
        echo "完成: ${sample}"
        ((done_n++))
    else
        echo "失败: ${sample}，请检查日志 ${sample_log}"
        ((fail_n++))
    fi
done < <(
    find . -maxdepth 1 -type f \
        \( -name '*_1.fq.gz' -o -name '*_1.fastq.gz' -o -name '*_1.clean.fq.gz' -o -name '*_1.sub10.fq.gz' \) \
        -print0 | sort -z
)

if [[ ${found} -eq 0 ]]; then
    echo "当前目录未找到任何 *_1.fq.gz / *_1.fastq.gz 类型的双端输入文件。"
    exit 1
fi

echo
echo "全部任务结束"
echo "成功: ${done_n}"
echo "跳过: ${skip_n}"
echo "失败: ${fail_n}"

if [[ ${fail_n} -gt 0 ]]; then
    exit 1
fi
