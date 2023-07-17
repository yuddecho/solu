# coding: utf-8

import os

is_test = False

root = '/public/home/yudong/protein'

if is_test:
    root = '/Users/yudd/code/python/echo/data/db'

log_file = f'{root}/train.log'

# train_solu_dataset = f'{root}/finally/tt_pdb_solu.fasta'
# train_insolu_dataset = f'{root}/finally/tt_insolu.fasta'
train_solu_dataset = f'{root}/finally/train_tt_pdb_solu.fasta'
train_insolu_dataset = f'{root}/finally/train_tt_insolu.fasta'

test_solu_dataset = f'{root}/finally/test_tt_pdb_solu.fasta'
test_insolu_dataset = f'{root}/finally/test_tt_insolu.fasta'

val_chang_solu = f'{root}/finally/chang_solu.fasta'
val_chang_insolu = f'{root}/finally/chang_insolu.fasta'

val_nesg_solu = f'{root}/finally/nesg_solu.fasta'
val_nesg_insolu = f'{root}/finally/nesg_insolu.fasta'

# test data info
# 34 1697
# 7906 7826.94 79.0600000000004
# 551
seq_max_len = 551
acid_numb = 34


if __name__ == '__main__':
    file_list = [train_solu_dataset, train_insolu_dataset]

    for file in file_list:
        if os.path.exists(file):
            print(file)

    print(f'Seq max len: {seq_max_len}')


