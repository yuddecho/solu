# coding: utf-8

import os

is_test = True

root = '/public/home/yudong/protein'

if is_test:
    root = '/Users/yudd/code/python/echo/data/db'

log_file = f'{root}/train.log'

train_solu_dataset = f'{root}/finally/tt_pdb_solu.fasta'
train_insolu_dataset = f'{root}/finally/tt_insolu.fasta'

test_chang_solu = f'{root}/finally/chang_solu.fasta'
test_chang_insolu = f'{root}/finally/chang_insolu.fasta'

test_nesg_solu = f'{root}/finally/nesg_solu.fasta'
test_nesg_insolu = f'{root}/finally/nesg_insolu.fasta'

seq_max_len = 1186
acid_numb = 21


if __name__ == '__main__':
    file_list = [train_solu_dataset, train_insolu_dataset]

    for file in file_list:
        if os.path.exists(file):
            print(file)

    print(f'Seq max len: {seq_max_len}')


