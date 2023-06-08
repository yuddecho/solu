# coding: utf-8

import os

is_test = False

root = '../data'

sample_train_solu_set = os.path.join(root, 'dataset/sample_train_solu_set.fasta')
sample_train_insolu_set = os.path.join(root, 'dataset/sample_train_insolu_set.fasta')

if is_test:
    train_solu_dataset = sample_train_solu_set
    train_insolu_dataset = sample_train_insolu_set
else:
    train_solu_dataset = f'{root}/dataset/solu.fasta'
    train_insolu_dataset = f'{root}/dataset/insolu.fasta'

seq_max_len = 1284


if __name__ == '__main__':
    file_list = [train_solu_dataset, train_insolu_dataset]

    for file in file_list:
        if os.path.exists(file):
            print(file)

    print(f'Seq max len: {seq_max_len}')


