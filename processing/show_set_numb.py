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


def get_tmbed(root_dir, file_names):
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


# cols = ['tt', 'tt_solu', 'tt_insolu', 'pdb_solu', 'chang', 'chang_solu', 'chang_insolu', 'nesg', 'nesg_solu', 'nesg_insolu']
#
# with open(f'{root}/db/res.csv', 'w', encoding='utf-8') as w:
#     w.write(f',,ttdb,tt-solu,tt-insolu,pdb,chang,chang-solu,chang-insolu,nesg,nesg-solu,nesg-insolu\n')
#     w.write(f',,ttdb,tt-solu,tt-insolu,pdb,chang,chang-solu,chang-insolu,nesg,nesg-solu,nesg-insolu\n')

# tmbed
files = ['tt_pdb_solu', 'tt_insolu', 'chang_solu', 'chang_insolu', 'nesg_solu', 'nesg_insolu']
tmbed_number = get_tmbed(f'{root}/db/tmbed', files)
print(tmbed_number)

# finally 裁剪长度
files = ['tt_pdb_solu', 'tt_insolu', 'chang_solu', 'chang_insolu', 'nesg_solu', 'nesg_insolu']
finally_number = {}
# seq_min, seq_max = 30, 622
seq_min, seq_max = 20, 1186
seq_count = 0

for fasta_file_name in files:
    fasta_file = f'{root}/db/tmbed/{fasta_file_name}.fasta'
    target_file = f'{root}/db/finally/{fasta_file_name}.fasta'

    with open(target_file, 'w', encoding='utf-8') as w:
        with open(fasta_file, 'r', encoding='utf-8') as r:
            while True:
                fasta = r.readline()
                seq = r.readline()

                seq = seq.strip()

                if not seq:
                    finally_number[fasta_file_name] = seq_count
                    seq_count = 0
                    break

                # 限制长度
                seq_len = len(seq)
                if seq_len < seq_min or seq_len > seq_max:
                    continue

                seq_count += 1

                w.write(fasta)
                w.write(f'{seq}\n')

print(finally_number)
