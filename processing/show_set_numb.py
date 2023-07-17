import os


# source
def get_source(root_dir):
    source_dir = f'{root_dir}/db/source'
    source_file = ['tt.csv', 'pdb_seqres.txt', 'chang_solu.csv', 'chang_insolu.csv', 'nesg.fasta']
    source_name = {'tt.csv': 'tt', 'pdb_seqres.txt': 'pdb_solu', 'chang_solu.csv': 'chang_solu', 'chang_insolu.csv': 'chang_insolu', 'nesg.fasta': 'nesg'}

    res = {}

    for file in source_file:
        file_path = f'{source_dir}/{file}'

        with open(file_path, 'r', encoding='utf-8') as r:
            if file == 'tt.csv':
                res[source_name[file]] = len(r.readlines())
            else:
                res[source_name[file]] = int(len(r.readlines()) / 2)

    return res


def get_cleaned(root_dir, cleaned_files):
    cleaned_dir = f'{root_dir}/db/cleaned'
    res = {}

    for file_name in cleaned_files:
        file_path = f'{cleaned_dir}/{file_name}.fasta'

        with open(file_path, 'r', encoding='utf-8') as r:
            res[file_name] = int(len(r.readlines()) / 2)

    return res


# mmseqs
def get_mmseqs(root_dir, mmseqs_files, tag):
    mmseqs_dir = f'{root_dir}/db/mmseqs'
    res = {}

    # 90
    for file_name in mmseqs_files:
        file_path = f'{mmseqs_dir}/{file_name}{tag}.fasta'

        if os.path.exists(file_path):
            with open(file_path, 'r', encoding='utf-8') as r:
                res[file_name] = int(len(r.readlines()) / 2)

    return res


def get_numbers(root_dir, file_names):
    res = {}

    for file_name in file_names:
        file_path = f'{root_dir}/{file_name}.fasta'

        with open(file_path, 'r', encoding='utf-8') as r:
            res[file_name] = int(len(r.readlines()) / 2)

    return res


root = '../../data'
files = ['tt_solu', 'tt_insolu', 'pdb_solu', 'chang_solu', 'chang_insolu', 'nesg_solu', 'nesg_insolu']

# source: ['tt.csv', 'pdb_seqres.txt', 'chang_solu.fasta', 'chang_insolu.fasta', 'nesg.fasta']
source_file_number = get_source(root)

# cleaned
cleaned_file_number = get_cleaned(root, files)

# mmseqs
"""
    各个数据源去除同一性超过90%的数据
    同一数据源可溶与不可溶数据之间去除同一性为 100% 的数据
    
    去除 tt_solu 中与 pdb_solu 同一性超过 90% 的数据
    
    去除 tt_insolu 中与 pdb_solu 同一性为 100% 的数据
    去除 chang_insolu 中与 pdb_solu 同一性为 100% 的数据
    去除 nesg_insolu 中与 pdb_solu 同一性为 100% 的数据
    
    合并 tt_solu 与 pdb_solu 数据，得到 tt_pdb_solu
    
    去除 tt_pdb_solu 中与 chang_solu 同一性超过 30% 的数据
    去除 tt_insolu 中与 chang_insolu 同一性超过 30% 的数据
    去除 tt_pdb_solu 中与 nesg_solu 同一性超过 30% 的数据
    去除 tt_insolu 中与 nesg_insolu 同一性超过 30% 的数据
"""
mmseqs_res_seq = get_mmseqs(root, files, "_rep_seq")
mmseqs_res_seq_100 = get_mmseqs(root, files, "_rep_seq_100")
mmseqs_res_seq_100_100 = get_mmseqs(root, files, "_rep_seq_100_100")
mmseqs_res_seq_100_100_30 = get_mmseqs(root, files, "_rep_seq_100_100_30")

print(source_file_number)
print(cleaned_file_number)

print(mmseqs_res_seq)
print(mmseqs_res_seq_100)
print(mmseqs_res_seq_100_100)
print(mmseqs_res_seq_100_100_30)

# tmbed
files = ['tt_pdb_solu', 'tt_insolu', 'chang_solu', 'chang_insolu', 'nesg_solu', 'nesg_insolu']
tmbed_number = get_numbers(f'{root}/db/tmbed', files)

# tt_pdb_solu 找出 tt 和 pdb
tt_solu_total, pdb_solu_total = 0, 0
with open(f'{root}/db/tmbed/tt_pdb_solu.fasta', 'r', encoding='utf-8') as r:
    while True:
        fasta = r.readline()
        seq = r.readline()

        fasta = fasta.strip()
        if not fasta:
            break

        if '_solu_pdb' in fasta:
            pdb_solu_total += 1
            continue

        if '_solu_tt' in fasta:
            tt_solu_total += 1

tmbed_number['tt_solu'], tmbed_number['pdb_solu'] = tt_solu_total, pdb_solu_total

print(tmbed_number)

# finally 裁剪长度之后
files = ['tt_pdb_solu', 'tt_insolu', 'chang_solu', 'chang_insolu', 'nesg_solu', 'nesg_insolu']
finally_number = get_numbers(f'{root}/db/finally', files)

# tt_pdb_solu 找出 tt 和 pdb
tt_solu_total, pdb_solu_total = 0, 0
with open(f'{root}/db/finally/tt_pdb_solu.fasta', 'r', encoding='utf-8') as r:
    while True:
        fasta = r.readline()
        seq = r.readline()

        fasta = fasta.strip()
        if not fasta:
            break

        if '_solu_pdb' in fasta:
            pdb_solu_total += 1
            continue

        if '_solu_tt' in fasta:
            tt_solu_total += 1

finally_number['tt_solu'], finally_number['pdb_solu'] = tt_solu_total, pdb_solu_total

print(finally_number)

cols = ['pdb_solu', 'tt', 'tt_solu', 'tt_insolu', 'chang', 'chang_solu', 'chang_insolu', 'nesg', 'nesg_solu', 'nesg_insolu']
data_dict = [source_file_number, cleaned_file_number, mmseqs_res_seq, mmseqs_res_seq_100, mmseqs_res_seq_100_100, mmseqs_res_seq_100_100_30, tmbed_number, finally_number]

with open(f'{root}/db/res.csv', 'w', encoding='utf-8') as w:
    for col_name in cols:
        w.write(f'{col_name},')
    w.write('\n')

    for data_item in data_dict:
        for col_name in cols:
            w.write(f'{data_item.get(col_name, "-")},')
        w.write('\n')























