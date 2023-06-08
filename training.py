# -*- coding: utf-8 -*-

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

import pickle
import os

from dataset import ProteinDataset
from models import KMersCNN
from args import root, train_solu_dataset, train_insolu_dataset


# 训练
def model_train(pds, te, bs, lr, mp, de):
    protein_dataset = pds
    train_epoch = te
    batch_size = bs
    model_pth_save_path = mp
    device = de

    # 数据
    train_loader = DataLoader(dataset=protein_dataset, batch_size=batch_size, shuffle=True)

    # 模型
    model = KMersCNN()
    if torch.cuda.device_count() > 1:
        print("Let's use", torch.cuda.device_count(), "GPUs!")
        # dim = 0 [30, xxx] -> [10, ...], [10, ...], [10, ...] on 3 GPUs
        model = nn.DataParallel(model)

    model.to(device)

    # 定义优化器和损失函数
    optimization = torch.optim.Adam(model.parameters(), lr)
    loss_func = nn.CrossEntropyLoss()  # 交叉熵

    loss_data = []
    before_loos = 10000

    # 训练网络
    for epoch in range(train_epoch):
        for step, (x, y) in enumerate(train_loader):
            x, y = x.to(device), y.to(device)

            out = model(x)
            loss = loss_func(out, y)

            optimization.zero_grad()
            loss.backward()
            optimization.step()

            # 每50步打印一下结果
            # if step % 2 == 0:
            loss_data.append(loss.item())
            print(f'Epoch:{epoch + 1} Step:{step + 1} Train loss:{loss.item()}')
            # Epoch:10 Step:600 Train loss:0.015473434701561928

            if loss.item() < before_loos:
                before_loos = loss.item()

                # 保存模型参数
                torch.save(model.state_dict(), model_pth_save_path)

    return loss_data


# main
def main():
    # 数据
    train_solu = train_solu_dataset
    train_insolu = train_insolu_dataset

    # args
    os.environ['CUDA_VISIBLE_DEVICES'] = '0,1,2,3'
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # acid_encoding = 'one-hot'
    acid_encoding = 'natural-number'
    train_epoch = 100
    batch_size = 200
    learn_rate = 0.01
    model_pth_save_path = f'{root}/model.pth'
    loss_pkl_save_path = f'{root}/loos.pkl'

    # dataset
    protein_dataset = ProteinDataset(train_solu, train_insolu, acid_encoding)

    # train
    loss = model_train(pds=protein_dataset, te=train_epoch, bs=batch_size, lr=learn_rate, mp=model_pth_save_path, de=device)

    # save loss
    with open(loss_pkl_save_path, 'wb') as w:
        pickle.dump(loss, w)


if __name__ == '__main__':
    main()
