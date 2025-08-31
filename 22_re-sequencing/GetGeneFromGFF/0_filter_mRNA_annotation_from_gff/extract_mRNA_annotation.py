def read_gff(filename):
    mrna_annotations = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # 忽略注释行
            fields = line.strip().split('\t')
            if fields[2] == 'mRNA':
                mrna_annotations.append(fields)
    return mrna_annotations

def remove_overlapping_annotations(annotations):
    non_overlapping = []
    sorted_annotations = sorted(annotations, key=lambda x: (x[0], int(x[3]), int(x[4])))
    current_chr, current_start, current_end = '', 0, 0
    for annotation in sorted_annotations:
        chr, start, end = annotation[0], int(annotation[3]), int(annotation[4])
        if chr != current_chr or start > current_end:
            non_overlapping.append(annotation)
            current_chr, current_start, current_end = chr, start, end
        else:
            current_end = max(current_end, end)
    return non_overlapping

def main():
    gff_file = 'ref.gff'
    mrna_annotations = read_gff(gff_file)
    filtered_annotations = remove_overlapping_annotations(mrna_annotations)

    # 输出结果
    for annotation in filtered_annotations:
        print('\t'.join(annotation))

if __name__ == '__main__':
    main()
