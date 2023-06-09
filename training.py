# -*- coding: utf-8 -*-

import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from torch.backends import cudnn

import pickle
import os
import time
import random
import numpy as np
from tqdm import tqdm

from dataset import ProteinDataset
from models import KMersCNN
from args import root, train_solu_dataset, train_insolu_dataset
from units import log


# main
class Training:
    def __init__(self):
        # args
        os.environ['CUDA_VISIBLE_DEVICES'] = '0,1,2,3'
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

        # 数据
        self.train_solu = train_solu_dataset
        self.train_insolu = train_insolu_dataset

        # acid_encoding = 'one-hot'
        self.acid_encoding = 'natural-number'
        self.epoch = 300
        self.curr_epoch = 1
        self.batch_size = 2000
        self.learn_rate = 0.001
        self.model_pth_save_path = f'{root}/model.pth'
        self.loss_pkl_save_path = f'{root}/loos.pkl'
        self.checkpoint_path = f'{root}/checkpoint.pt'

        self.model = KMersCNN()

        # 定义优化器和损失函数
        self.optimizer = torch.optim.Adam(self.model.parameters(), self.learn_rate)
        self.loss_func = nn.CrossEntropyLoss()  # 交叉熵

        self.loss_data = f'{root}/loos.scv'
        self.loss_data_a = open(self.loss_data, 'a', encoding='utf-8')
        self.last_loos = 10000

        # checkpoint
        self.resume = True
        if self.resume and os.path.exists(self.checkpoint_path):
            checkpoint = torch.load(self.checkpoint_path)
            self.model.load_state_dict(checkpoint['model_state_dict'])
            self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
            self.curr_epoch = checkpoint['epoch'] + 1
            self.last_loss = checkpoint['loss']

        self.last_stop_round_count_epoch = 3

    def print(self):
        # print log
        log(f'\n##############  dataset  ##############')
        log(f'solu: {self.train_solu}')
        log(f'insolu: {self.train_insolu}')
        log(f'encoding: {self.acid_encoding}')

        log(f'\n##############  model  ##############')
        log(f'epoch: {self.epoch}')
        log(f'batch size: {self.batch_size}')
        log(f'learn rate: {self.learn_rate}')
        log(f'device: {self.device}')
        log(f'{self.model_pth_save_path}')
        log(f'{self.loss_pkl_save_path}')
        log(f'\n')

    def run(self):
        # dataset
        protein_dataset = ProteinDataset(self.train_solu, self.train_insolu, self.acid_encoding)

        # dataloader
        train_loader = DataLoader(dataset=protein_dataset, batch_size=self.batch_size, shuffle=True)

        # 模型
        if torch.cuda.device_count() > 1:
            log(f"Let's use, {torch.cuda.device_count()}, GPUs!")
            # dim = 0 [30, xxx] -> [10, ...], [10, ...], [10, ...] on 3 GPUs
            self.model = nn.DataParallel(self.model)

        self.model.to(self.device)

        # 训练网络
        bar = tqdm(total=self.epoch, desc="started training")
        bar.update(self.curr_epoch)
        for epoch in range(self.curr_epoch, self.epoch):
            loss_data = []
            for step, (x, y) in enumerate(train_loader):
                x, y = x.to(self.device), y.to(self.device)

                out = self.model(x)
                loss = self.loss_func(out, y)

                self.optimizer.zero_grad()
                loss.backward()
                self.optimizer.step()

                # 每50步打印一下结果
                # if step % 2 == 0:
                loss_data.append(loss.item())
                self.loss_data_a.write(f'{epoch},{step},{loss.item()}\n')
                log(f'Epoch:{epoch + 1} Step:{step + 1} Train loss:{loss.item()}', is_print=False)

                if loss.item() < self.last_loos:
                    self.last_stop_round_count_epoch = 3
                    self.last_loos = loss.item()

                    # 保存模型参数
                    checkpoint = {
                        'epoch': epoch,
                        'model_state_dict': self.model.state_dict(),
                        'optimizer_state_dict': self.optimizer.state_dict(),
                        'loss': self.last_loos,
                    }
                    torch.save(checkpoint, self.model_pth_save_path)

            bar.update(1)
            bar.set_description(f'epoch: {epoch}, loos: {sum(loss_data)/len(loss_data)}')

            self.last_stop_round_count_epoch -= 1
            if self.last_stop_round_count_epoch == 0:
                break

    def __del__(self):
        self.loss_data_a.close()


def setup_seed(seed=2023):
    random.seed(seed)  # Python的随机性
    np.random.seed(seed)  # numpy的随机性
    os.environ['PYTHONHASHSEED'] = str(seed)  # 设置Python哈希种子，为了禁止hash随机化，使得实验可复现

    torch.manual_seed(seed)  # torch的CPU随机性，为CPU设置随机种子
    torch.cuda.manual_seed(seed)  # torch的GPU随机性，为当前GPU设置随机种子
    torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.   torch的GPU随机性，为所有GPU设置随机种子

    # 这俩严重降低效率
    # torch.backends.cudnn.deterministic = True  # 选择确定性算法
    # torch.backends.cudnn.benchmark = False  # if benchmark=True, deterministic will be False

    torch.backends.cudnn.enabled = False

    # 尽管设置的了种子，但nn.Embedding()初始化还是会不一致，可使用持久化保持


if __name__ == '__main__':
    setup_seed()

    training = Training()
    training.print()
    training.run()
