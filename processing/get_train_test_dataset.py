import random

is_test = False


def dataset_random(source_file, train_target_file, test_target_file, source_num, target_num):
    index_seq = list(range(int(source_num)))
    random.shuffle(index_seq)

    target_seq = index_seq[:target_num]

    print(target_seq)

    index_count = 0
    with open(train_target_file, 'w', encoding='utf-8') as w_train:
        with open(test_target_file, 'w', encoding='utf-8') as w_test:
            with open(source_file, 'r', encoding='utf-8') as r:
                while True:
                    fasta = r.readline()
                    seq = r.readline()

                    seq = seq.strip()
                    if not seq:
                        break

                    if index_count in target_seq:
                        w_test.write(fasta)
                        w_test.write(f'{seq}\n')
                    else:
                        w_train.write(fasta)
                        w_train.write(f'{seq}\n')

                    index_count += 1


def get_line_num(file_name):
    count = 0
    with open(file_name, 'r', encoding='utf-8') as r:
        rows = r.readlines()
        for row in rows:
            row = row.strip()
            if not row:
                break

            count += 1

    return count


if is_test:
    root = '../../data/db'
else:
    root = '/public/home/yudong/protein'

solu_dataset = f'{root}/finally/tt_pdb_solu.fasta'
insolu_dataset = f'{root}/finally/tt_insolu.fasta'

train_solu_dataset = f'{root}/finally/train_tt_pdb_solu.fasta'
train_insolu_dataset = f'{root}/finally/train_tt_insolu.fasta'

test_solu_dataset = f'{root}/finally/test_tt_pdb_solu.fasta'
test_insolu_dataset = f'{root}/finally/test_tt_insolu.fasta'

# 划分数据集
test_num = 5000

solu_num = get_line_num(solu_dataset) / 2
insolu_num = get_line_num(insolu_dataset) / 2

print(solu_num, insolu_num)

dataset_random(solu_dataset, train_solu_dataset, test_solu_dataset, solu_num, test_num)
dataset_random(insolu_dataset, train_insolu_dataset, test_insolu_dataset, insolu_num, test_num)


