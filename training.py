# -*- coding: utf-8 -*-

import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from torch.backends import cudnn

import os
import random
import numpy as np
from tqdm import tqdm

from dataset import ProteinDataset
from models import KMersCNN
from args import root, train_solu_dataset, train_insolu_dataset, test_chang_solu, test_chang_insolu, test_nesg_solu, test_nesg_insolu, is_test
from units import log


# main
class Training:
    def __init__(self, resume):
        # args
        os.environ['CUDA_VISIBLE_DEVICES'] = '0,1,2,3'
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

        # 数据
        self.train_solu = train_solu_dataset
        self.train_insolu = train_insolu_dataset

        self.test_chang_solu = test_chang_solu
        self.test_chang_insolu = test_chang_insolu

        self.test_nesg_solu = test_nesg_solu
        self.test_nesg_insolu = test_nesg_insolu

        # acid_encoding = 'one-hot'
        self.acid_encoding = 'natural-number'
        self.epoch = 900
        self.curr_epoch = 0
        self.batch_size = 1024
        if is_test:
            self.batch_size = 2
        self.learn_rate = 0.001

        self.checkpoint_path = f'{root}/checkpoint.pt'

        self.model = KMersCNN()

        # 定义优化器和损失函数
        self.optimizer = torch.optim.Adam(self.model.parameters(), self.learn_rate)
        self.loss_func = nn.CrossEntropyLoss()  # 交叉熵

        self.loss_data = f'{root}/loos.csv'

        self.last_acc = -1

        # checkpoint
        self.resume = resume
        if self.resume and os.path.exists(self.checkpoint_path):
            checkpoint = torch.load(self.checkpoint_path)
            self.model.load_state_dict(checkpoint['model_state_dict'])
            self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
            self.curr_epoch = checkpoint['epoch'] + 1
            self.last_acc = checkpoint['acc']

    def print(self):
        # print log
        log(f'##############  dataset  ##############')
        log(f'solu: {self.train_solu}')
        log(f'insolu: {self.train_insolu}')
        log(f'encoding: {self.acid_encoding}')

        log(f'\n##############  model  ##############')
        log(f'epoch: {self.epoch}')
        log(f'batch size: {self.batch_size}')
        log(f'learn rate: {self.learn_rate}')
        log(f'device: {self.device}')
        log(f'{self.checkpoint_path}')
        log(f'{self.loss_data}')
        log(f'\n')

    def save_loss(self, msg):
        with open(self.loss_data, 'a', encoding='utf-8') as w:
            w.write(f'{msg}\n')

    def _train(self, epoch, data_loader):
        # https://www.jianshu.com/p/1cc53eef82bf
        # 每一批训练
        self.model.train()
        total_loss_data = 0
        total_acc = 0

        with torch.enable_grad():
            for step, (x, y) in enumerate(data_loader):
                x, y = x.to(self.device), y.to(self.device)

                out = self.model(x)
                loss = self.loss_func(out, y)

                acc = ((out.argmax(1) == y).sum().item()) / len(y)

                loss_value = loss.item()

                total_acc += acc
                total_loss_data += loss_value

                self.optimizer.zero_grad()
                loss.backward()
                self.optimizer.step()

                log(f'Info: train Epoch {epoch}, Step {step}, Acc {acc}, Loss {loss_value}', False)

        return total_loss_data / len(data_loader), total_acc / len(data_loader)

    def _test(self, epoch, data_loader, test_set_name):
        self.model.eval()
        total_loss_data = 0
        total_acc = 0

        with torch.no_grad():
            for step, (x, y) in enumerate(data_loader):
                x, y = x.to(self.device), y.to(self.device)

                out = self.model(x)
                loss = self.loss_func(out, y)

                acc = ((out.argmax(1) == y).sum().item()) / len(y)
                loss_value = loss.item()

                total_acc += acc
                total_loss_data += loss_value

                log(f'Info: Test {test_set_name} Epoch {epoch}, Step {step}, Acc {acc}, Loss {loss_value}', False)

        return total_loss_data / len(data_loader), total_acc / len(data_loader)

        pass

    def run(self):
        # dataset
        protein_dataset = ProteinDataset(self.train_solu, self.train_insolu, self.acid_encoding)
        chang_dataset = ProteinDataset(self.test_chang_solu, self.test_chang_insolu, self.acid_encoding)
        nesg_dataset = ProteinDataset(self.test_nesg_solu, self.test_nesg_insolu, self.acid_encoding)

        # dataloader
        train_loader = DataLoader(dataset=protein_dataset, batch_size=self.batch_size, shuffle=True)
        chang_loader = DataLoader(dataset=chang_dataset, batch_size=self.batch_size, shuffle=True)
        nesg_loader = DataLoader(dataset=nesg_dataset, batch_size=self.batch_size, shuffle=True)

        # 模型
        if torch.cuda.device_count() > 1:
            log(f"Let's use, {torch.cuda.device_count()}, GPUs!")
            # dim = 0 [30, xxx] -> [10, ...], [10, ...], [10, ...] on 3 GPUs
            self.model = nn.DataParallel(self.model)

        self.model.to(self.device)

        # 训练网络
        bar = tqdm(total=self.epoch-1, desc="started training")
        bar.update(self.curr_epoch)
        for epoch in range(self.curr_epoch, self.epoch):
            # train
            loss_train_data, acc_train = self._train(epoch, train_loader)

            # test
            loss_chang_data, chang_acc = self._test(epoch, chang_loader, 'chang')
            loss_nesg_data, nesg_acc = self._test(epoch, nesg_loader, 'nesg')

            self.save_loss(f'{loss_train_data},{loss_chang_data},{loss_nesg_data},{acc_train},{chang_acc},{nesg_acc}')

            acc = (nesg_acc + chang_acc) / 2

            # 保存数据
            if acc > self.last_acc:
                log(f'Info: save model, last-acc {self.last_acc}, curr-acc {acc} > {loss_train_data},{loss_chang_data},{loss_nesg_data},{acc_train},{chang_acc},{nesg_acc}')

                self.last_acc = acc

                # 保存模型参数
                checkpoint = {
                    'epoch': epoch,
                    'model_state_dict': self.model.state_dict(),
                    'optimizer_state_dict': self.optimizer.state_dict(),
                    'loss': self.last_acc,
                }

                torch.save(checkpoint, self.checkpoint_path)

            bar.update(1)
            bar.set_description(f'epoch: {epoch}, acc: {acc} loss: {(loss_nesg_data + loss_chang_data) / 2}')


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


def init(resume):
    setup_seed()

    # 删除文件
    if not resume:
        loos_file = f'{root}/loos.csv'
        if os.path.exists(loos_file):
            os.remove(loos_file)

        log_file = f'{root}/train.log'
        if os.path.exists(log_file):
            os.remove(log_file)


if __name__ == '__main__':
    is_resume = False
    init(is_resume)

    training = Training(is_resume)
    training.print()
    training.run()
