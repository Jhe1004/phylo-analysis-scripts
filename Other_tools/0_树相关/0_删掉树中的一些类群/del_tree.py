import os
import shutil
from ete3 import Tree

# --- 配置 ---
now_dir = os.getcwd()
del_file_list_file = "tree.txt"
output_dir_name = "processed_trees"  # 输出目录名称
output_dir_path = os.path.join(now_dir, output_dir_name) # 输出目录完整路径
file_prefix = "dating." # 定义文件名前缀

# --- 函数 ---
def get_tree_files(directory, prefix):
    """查找指定目录中所有以指定前缀开头的文件。"""
    file_list = []
    try:
        for each in os.listdir(directory):
            # 确保是文件，并且以指定前缀开头
            if os.path.isfile(os.path.join(directory, each)) and each.startswith(prefix):
                file_list.append(each)
    except FileNotFoundError:
        print(f"错误：找不到目录 '{directory}'")
    return file_list

# --- 主程序逻辑 ---
if __name__ == "__main__":
    # 1. 创建输出目录 (如果不存在)
    try:
        os.makedirs(output_dir_path, exist_ok=True)
        print(f"输出目录: '{output_dir_path}' 已确保存在。")
    except OSError as e:
        print(f"错误：无法创建目录 '{output_dir_path}': {e}")
        exit()

    # 2. 获取当前目录下符合前缀的树文件列表
    tree_files_in_current_dir = get_tree_files(now_dir, file_prefix) # 使用定义的前缀

    if not tree_files_in_current_dir:
        # 更新提示信息，说明查找的是特定前缀的文件
        print(f"在 '{now_dir}' 中未找到以 '{file_prefix}' 开头的文件。")
        exit()

    # 3. 读取需要删除的叶子节点列表 (只读一次)
    del_leaves_set = set()
    try:
        with open(del_file_list_file, "r") as read_file:
            for each_line in read_file:
                del_leaves_set.add(each_line.strip())
        print(f"从 '{del_file_list_file}' 读取了 {len(del_leaves_set)} 个要删除的叶子节点名称。")
    except FileNotFoundError:
        print(f"警告：删除列表文件 '{del_file_list_file}' 未找到。将不会删除任何叶子节点。")
        # If the delete list is essential, you might want to exit() here
        # exit()


    # 4. 将所有找到的原始树文件复制到输出目录
    print(f"\n正在将 {len(tree_files_in_current_dir)} 个以 '{file_prefix}' 开头的原始文件复制到 '{output_dir_name}'...")
    copy_errors = 0
    for tree_file in tree_files_in_current_dir:
        source_path = os.path.join(now_dir, tree_file)
        destination_path = os.path.join(output_dir_path, tree_file)
        try:
            shutil.copy2(source_path, destination_path) # copy2 preserves metadata like modification time
        except Exception as e:
            print(f"  复制文件 '{tree_file}' 时出错: {e}")
            copy_errors += 1
    if copy_errors == 0:
        print("所有文件复制完成。")
    else:
        print(f"{copy_errors} 个文件复制失败。")


    # 5. 处理每个树文件，如有必要则在输出目录中覆盖更新
    print(f"\n正在处理 {len(tree_files_in_current_dir)} 个树文件...")
    processed_count = 0
    skipped_count = 0
    error_count = 0

    for tree_file in tree_files_in_current_dir:
        original_file_path = os.path.join(now_dir, tree_file) # Path to original file (for loading)
        output_file_path = os.path.join(output_dir_path, tree_file) # Path to file in output dir (for writing/overwriting)

        try:
            # 加载原始树文件 (RAxML bipartition 文件通常是 Newick 格式)
            # format=1 通常适用于标准 Newick，但如果遇到问题可能需要尝试其他格式 (0-9)
            # Ensure you load the original file, not the potentially already copied one if reprocessing
            t = Tree(original_file_path, format=1)

            original_leaves = set(t.get_leaf_names())
            if not original_leaves:
                print(f"  跳过 {tree_file}: 未找到叶子节点。")
                skipped_count += 1
                continue

            # 确定要保留的叶子节点
            leaves_to_keep = [name for name in original_leaves if name not in del_leaves_set]

            # 检查是否有叶子节点被实际移除
            if len(leaves_to_keep) < len(original_leaves):
                num_removed = len(original_leaves) - len(leaves_to_keep)
                print(f"  处理中 {tree_file}: 准备移除 {num_removed} 个叶子节点...")

                if not leaves_to_keep:
                    print(f"    警告：文件 {tree_file} 中的所有叶子节点都将被删除。将生成一个空树或只包含根节点的树。")
                    # Decide how to handle this - write an empty file, skip, or let prune handle it?
                    # ETE's prune might raise an error if leaves_to_keep is empty.
                    # Let's add a check:
                    if not leaves_to_keep:
                         print(f"    跳过写入空树: {tree_file}")
                         # Optionally delete the copied file if you don't want empty results
                         # try:
                         #     os.remove(output_file_path)
                         # except OSError:
                         #     pass
                         skipped_count += 1 # Count as skipped if we don't write
                         continue # Skip to next file


                # --- POTENTIAL ISSUE FOR ULTRAMETRIC TREES ---
                # The prune function recalculates branch lengths (`dist`) based on the
                # remaining topology. This preserves relative evolutionary distances
                # but often breaks the fixed node ages (heights) required for
                # ultrametricity. There's no simple flag in ete3's prune
                # to force preservation of original node ages.
                t.prune(leaves_to_keep)
                # --- END POTENTIAL ISSUE ---

                # 将修剪后的树写回到输出目录中的同名文件 (覆盖复制的那个)
                # Use format=1 again, assuming Newick output is desired.
                # Using the standard format (often 0 or 1) usually preserves branch lengths
                # as calculated by prune.
                t.write(format=1, outfile=output_file_path)
                print(f"  已处理并更新: {tree_file} 在 '{output_dir_name}' 目录中。")
                processed_count += 1
            else:
                # If no leaves needed removal, the originally copied file is correct.
                print(f"  跳过修剪 {tree_file}: 无需删除叶子节点 (原始文件已复制到 '{output_dir_name}')。")
                skipped_count += 1

        except Exception as e:
            print(f"  处理文件 '{tree_file}' 时出错: {e}")
            # Consider logging the specific error for debugging
            # import traceback
            # print(traceback.format_exc())
            error_count += 1

    # --- 结束总结 ---
    print("\n--- 处理完成 ---")
    print(f"总共检查以 '{file_prefix}' 开头的文件数: {len(tree_files_in_current_dir)}") # 更新总结信息
    print(f"已修剪并更新的文件数: {processed_count}")
    print(f"无需修剪的文件数 (已复制): {skipped_count}")
    print(f"处理过程中出错的文件数: {error_count}")
    print(f"所有结果文件（无论是否修剪）均位于: '{output_dir_path}'")
    if error_count > 0:
        print(f"警告：{error_count} 个文件在处理过程中遇到错误。请检查上面的日志。")
    if processed_count > 0 :
         print("注意：如果输入树是超度量树（时间树），标准的剪枝操作可能已破坏其超度量属性。")