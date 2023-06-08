# coding: utf-8

import matplotlib.pyplot as plt


# 数据分布曲线，用来查看loss的变化
def polt_curve(data):
    # 数据分布曲线
    fig = plt.figure()
    plt.plot(range(len(data)), data, color='blue')
    plt.legend(['value'], loc='upper right')
    plt.xlabel('step')
    plt.ylabel('value')
    plt.show()


# 二分查找本地one-hot已编码文件
def binary_search(strs: list, target: str) -> bool:
    # one-hot 编码时，有超过85G数据，所以将编码的内容存为Python对象文件，该函数用于查询是否已编码，减少重复工作
    left, middle = 0, 0
    right = len(strs)

    if right == 0:
        print(f'Error: strs len is 0')
        raise

    while left < right:
        middle = left + ((right - left) >> 1)

        if strs[middle] > target:
            right = middle
            continue

        if strs[middle] < target:
            left = middle + 1
            continue

        if strs[middle] == target:
            return True

    return False
