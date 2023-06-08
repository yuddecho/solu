#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os  # os为内置模块，是一定存在的

while True:
    try:
        import torch
        from tqdm import tqdm
        import numpy
        import matplotlib.pyplot as plt

        print('All installation packages have been installed')
        break

    except ImportError as e:
        print(f'{e.msg}, Preparing to start installing {e.name}...')

        os.system(f'conda install {e.name} -y')

