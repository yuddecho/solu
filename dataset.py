# -*- coding: utf-8 -*-

from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

import numpy as np
import pickle
import os

from units import binary_search, log
from args import root, train_solu_dataset, train_insolu_dataset, seq_max_len, acid_numb, is_test


# ProteinDataset class
class ProteinDataset(Dataset):
    """ 蛋白质数据集 """

    def __init__(self, protein_solu_fasta, protein_insolu_fasta,
                 acid_encoding):
        """
        protein_solu_fasta: 溶解数据 fasta 格式文件
        protein_insolu_fasta: 不溶数据 fasta 格式文件
        """

        # 从根本区分
        # name, sequence, label, seq-encode
        self.id = []
        self.id_seq = {}
        self.id_label = {}

        # 序列数目
        # self.seq_num = 20000

        # 序列最大长度
        self.seq_max_len = seq_max_len

        # 氨基酸种类
        self.acid_class_num = acid_numb

        # 氨基酸编码方式
        '''
        acid_encoding:
            'one-hot',
            'natural-number'
        '''
        self.acid_encoding = acid_encoding

        # id: natural-number
        self.id_seq_natural_number_encode = {}
        self.seq_natural_number_recode = False  # False 尽可能是读取已编码数据
        self.id_seq_natural_number_encode_dict_pkl = f'{root}/seq_nn_encode.pkl'
        # log(f'{self.id_seq_natural_number_encode_dict_pkl}')

        # one hot ecode file for dir : dict: python pickle object, key: id, value: seq encode vector
        self.seq_one_hot_encode_dir = f'{root}/one_hot'
        # log(f'{self.seq_one_hot_encode_dir}')

        # init

        # read data fun
        self._read_fasta_file(protein_solu_fasta, protein_insolu_fasta)

        # seq encode fun
        self._sequence_encode()

    def _read_fasta_file(self, protein_solu_fasta, protein_insolu_fasta):
        """
        读取 fasta 文件
        """
        test_seq_num = 0
        for fasta_file in [protein_solu_fasta, protein_insolu_fasta]:
            with open(fasta_file, 'r', encoding='utf-8') as r:
                while True:
                    fasta = r.readline()
                    seq = r.readline()
                    seq = seq.strip()

                    if not seq:
                        break

                    # label
                    fasta = fasta.strip()
                    fasta = fasta[1:]

                    label = -1
                    if '_solu_' in fasta:
                        label = 1
                    if '_insolu_' in fasta:
                        label = 0

                    if label != -1:
                        self.id.append(f'{fasta}')
                        self.id_seq[fasta] = seq
                        self.id_label[fasta] = label

                        # 测试时，每个文件取 1000
                        test_seq_num += 1
                        if is_test and test_seq_num == 1000:
                            test_seq_num = 0
                            break

                    else:
                        print(fasta)

        self.seq_num = len(self.id)

    def _sequence_encode(self):
        if self.acid_encoding == 'one-hot':
            self._sequence_one_hot_encode()
            return

        if self.acid_encoding == 'natural-number':
            self._sequence_natural_number_encode()
            return

        print(f'Error: {self.acid_encoding} is not exists.')
        raise

    def _sequence_one_hot_encode(self):
        # one-hot
        acid_to_number = {
            'U': 0,
            'L': 1,
            'S': 2,
            'R': 3,
            'V': 4,
            'G': 5,
            'A': 6,
            'X': 7,
            'Q': 8,
            'I': 9,
            'O': 10,
            'D': 11,
            'P': 12,
            'E': 13,
            'W': 14,
            'H': 15,
            'T': 16,
            'Y': 17,
            'C': 18,
            'Z': 19,
            'K': 20,
            'M': 21,
            'N': 22,
            'F': 23,
            'B': 24
        }

        # 已编码文件合集
        existing_seq_encode_files = os.listdir(self.seq_one_hot_encode_dir)

        bar = tqdm(total=self.seq_num, desc='Sequence encode(oh)')
        for seq_id in self.id:
            # one_hot_encode_pkl = os.path.join(self.seq_one_hot_encode_dir, f'{seq_id}.pkl')
            one_hot_encode_pkl = f'{seq_id}.pkl'

            # 采用二分查找，从总体 8 分钟，降到 7 秒
            if binary_search(existing_seq_encode_files, one_hot_encode_pkl):
                bar.update(1)
                continue

            sequence = self.id_seq[seq_id]

            # 每条序列对应一个二维 编码矩阵 2004*25: seq_max_len * acid_class_num
            seq_item_encode = np.zeros([self.seq_max_len, self.acid_class_num])

            for seq_char_index, seq_char in enumerate(sequence):
                # 获取氨基酸对应位置下标
                seq_item_encode_acid_index = acid_to_number.get(seq_char, -1)

                # 对应位置赋值为1
                seq_item_encode[seq_char_index][seq_item_encode_acid_index] = 1

            # 完成一条序列编码, 写入文件
            # 一个文件 401KB, 总数据量预估 82GB
            with open(one_hot_encode_pkl, 'wb') as f:
                pickle.dump(seq_item_encode, f)

            bar.update(1)

        bar.close()

    def _sequence_natural_number_encode(self):
        # 读取已编码数据 215874条数据，3.4G大小
        # if not self.seq_natural_number_recode and os.path.exists(
        #         self.id_seq_natural_number_encode_dict_pkl):
        #     log(
        #         f'Info: natural number encode, find {self.id_seq_natural_number_encode_dict_pkl}, loading...'
        #     )
        #     with open(self.id_seq_natural_number_encode_dict_pkl, 'rb') as f:
        #         self.id_seq_natural_number_encode = pickle.load(f)
        #     return

        # natural number id
        acid_to_number = {
            'Y': 21,
            'H': 1,
            'F': 2,
            'M': 3,
            'W': 4,
            'S': 5,
            'R': 6,
            'L': 7,
            'A': 8,
            'G': 9,
            'I': 10,
            'N': 11,
            'K': 12,
            'D': 13,
            'V': 14,
            'E': 15,
            'X': 16,
            'C': 17,
            'Q': 18,
            'T': 19,
            'P': 20
        }

        bar = tqdm(total=self.seq_num, desc='Sequence encode(nn)')
        for seq_id in self.id:
            sequence = self.id_seq[seq_id]

            # 序列编码, 固定长度, 缺失补0，补的0会影响结果吗？
            # 每条序列对应一个一维编码矩阵 2004: seq_max_len
            seq_item_encode = np.zeros(self.seq_max_len)

            for seq_char_index, seq_char in enumerate(sequence):

                # 获取氨基酸对应位置下标, 即对应数字编码
                seq_item_encode_acid_index = acid_to_number.get(seq_char, -1)

                if seq_item_encode_acid_index == -1:
                    print(f'Error: {seq_id}, "{seq_char}" not define.')
                    raise

                # 对应位置赋值为1
                seq_item_encode[seq_char_index] = seq_item_encode_acid_index

            # 完成一条序列编码
            self.id_seq_natural_number_encode[seq_id] = seq_item_encode

            bar.update(1)

        # 将整体结果写入文件
        # with open(self.id_seq_natural_number_encode_dict_pkl, 'wb') as f:
        #     pickle.dump(self.id_seq_natural_number_encode, f)

        bar.close()

    def __getitem__(self, index):
        seq_id = self.id[index]
        encode_vector = None

        if self.acid_encoding == 'one-hot':
            one_hot_encode_pkl = os.path.join(self.seq_one_hot_encode_dir,
                                              f'{seq_id}.pkl')
            with open(one_hot_encode_pkl, 'rb') as f:
                encode_vector = pickle.load(f)

        if self.acid_encoding == 'natural-number':
            encode_vector = self.id_seq_natural_number_encode.get(seq_id, None)

        if encode_vector is None:
            raise

        return encode_vector, self.id_label[seq_id]

    def __len__(self):
        return self.seq_num


# Test ProteinDataset class
if __name__ == '__main__':
    train_solu = train_solu_dataset
    train_insolu = train_insolu_dataset

    # acid_encoding = 'one-hot'
    acids_encoding = 'natural-number'

    protein_dataset = ProteinDataset(train_solu, train_insolu, acids_encoding)
    print(len(protein_dataset))

    for i in range(2):
        print(protein_dataset[i])

    # 数据加载，每批次 4 个
    train_loader = DataLoader(dataset=protein_dataset, batch_size=4, shuffle=True)

    #
    for step, (b_x, b_y) in enumerate(train_loader):
        print(step)
        print(b_x)
        print(b_y)
        break
