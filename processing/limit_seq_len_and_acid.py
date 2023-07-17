# finally 裁剪长度
files = ['tt_pdb_solu', 'tt_insolu', 'chang_solu', 'chang_insolu', 'nesg_solu', 'nesg_insolu']
finally_number = {}
# seq_min, seq_max = 30, 622
seq_min, seq_max = 34, 551
seq_count = 0

is_test = False

if is_test:
    root = '../../data/db'
else:
    root = '/public/home/yudong/protein'

test_acid_set = ['H', 'F', 'M', 'W', 'S', 'R', 'L', 'A', 'G', 'I', 'N', 'K', 'D', 'V', 'E', 'X', 'C', 'Q', 'T', 'P', 'Y']
test_acid_set = set(test_acid_set)

for fasta_file_name in files:
    fasta_file = f'{root}/tmbed/{fasta_file_name}.fasta'
    target_file = f'{root}/finally/{fasta_file_name}.fasta'

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

                # 对训练集，限制氨基酸种类
                if fasta_file_name == 'tt_pdb_solu' or fasta_file_name == 'tt_insolu':
                    seq_set = set(seq)
                    # 差集
                    res_set = seq_set - test_acid_set
                    if res_set:
                        # 不为空
                        continue

                seq_count += 1

                w.write(fasta)
                w.write(f'{seq}\n')


# {'tt_pdb_solu': 49556, 'tt_insolu': 80084, 'chang_solu': 528, 'chang_insolu': 340, 'nesg_solu': 4356, 'nesg_insolu': 2487}
print(finally_number)
