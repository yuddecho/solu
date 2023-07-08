"""
1. 确定长度
1.5 氨基酸种类
2. 数据集平衡
"""
import os.path

import matplotlib.pyplot as plt
import pickle

root = '../../data'

tmbed_dir = f'{root}/db/tmbed'

train_set = ['tt_pdb_solu', 'tt_insolu']
test_set = ['chang_solu', 'chang_insolu', 'nesg_solu', 'nesg_insolu']

seqLen_seqNum_dict_pkl = f'{tmbed_dir}/seqLen_seqNum_dict.pkl'

seqLen_seqNum: dict = {}
acids_set = set()
seq_count = 0
re_do = True

if not os.path.exists(seqLen_seqNum_dict_pkl) or re_do:
    # 获取 测试集 长度分度
    for fasta_file in test_set:
        fasta_file = f'{tmbed_dir}/{fasta_file}.fasta'
        with open(fasta_file, 'r', encoding='utf-8') as r:
            while True:
                fasta = r.readline()
                seq = r.readline()

                seq = seq.strip()

                if not seq:
                    print(f'{fasta_file}: {seq_count}')
                    seq_count = 0
                    break

                # 统计氨基酸种类
                acids = set(seq)
                # 集合合并
                acids_set = acids_set.union(acids)

                # 统计序列长度
                seq_len = len(seq)
                seqLen_seqNum[seq_len] = seqLen_seqNum.get(seq_len, 0) + 1

                seq_count += 1

    with open(seqLen_seqNum_dict_pkl, 'wb') as f:
        pickle.dump(seqLen_seqNum, f)
else:
    with open(seqLen_seqNum_dict_pkl, 'rb') as f:
        seqLen_seqNum = pickle.load(f)

# 长度-数量 分布图 ==> 先确定为 20-1500, 测试集集中在 34-622
# 做一个分界线，分别占据多少数据量
# y 轴 1-100 放大，


seq_len_list = list(seqLen_seqNum.keys())

seq_len_list.sort()

seq_max_len = max(seq_len_list)
print(min(seq_len_list), seq_max_len)

x, y = [], []
seq_total = 0
for item in range(seq_max_len + 1):
    num_val = seqLen_seqNum.get(item, 0)
    if num_val != 0:
        x.append(item)
        y.append(num_val)
        seq_total += num_val

seq_total_present = seq_total * 0.99
print(seq_total, seq_total_present, seq_total - seq_total_present)

y_total = 0
x_limit_index = 0
for index, item in enumerate(y):
    y_total += item
    if y_total > seq_total_present:
        x_limit_index = index
        break
print(x_limit_index)

# y = [seqLen_seqNum.get(seq_len_item, 0) for seq_len_item in x]

plt.figure(1)
# plt.bar(x, y, facecolor='#9999ff', edgecolor='white')
plt.scatter(x, y, s=0.2, color='#9999ff')
# plt.plot(x, y)

# 画一条虚线
plt.plot([0, seq_max_len], [10, 10], 'k--')
# plt.plot([0, seq_max_len], [0, 0], 'k--')
plt.plot([x_limit_index, x_limit_index], [0, max(seqLen_seqNum.values())], 'k--')
# plt.plot([20, 20], [0, max(seqLen_seqNum.values())], 'k--')

# plt.show()
# plt.draw()

# 25 {'V', 'A', 'I', 'G', 'U', 'M', 'K', 'X', 'P', 'O', 'N', 'D', 'C', 'Y', 'Q', 'F', 'R', 'S', 'Z', 'L', 'E', 'W', 'T', 'B', 'H'}
# 21 {'H', 'F', 'M', 'W', 'S', 'R', 'L', 'A', 'G', 'I', 'N', 'K', 'D', 'V', 'E', 'X', 'C', 'Q', 'T', 'P', 'Y'}

print(len(acids_set), acids_set)

t1 = ['V', 'A', 'I', 'G', 'U', 'M', 'K', 'X', 'P', 'O', 'N', 'D', 'C', 'Y', 'Q', 'F', 'R', 'S', 'Z', 'L', 'E', 'W', 'T', 'B', 'H']
t2 = ['H', 'F', 'M', 'W', 'S', 'R', 'L', 'A', 'G', 'I', 'N', 'K', 'D', 'V', 'E', 'X', 'C', 'Q', 'T', 'P', 'Y']

set1 = set(t1)
set2 = set(t2)

t3 = set1 - set2
print(t3, len(t3))

for i in range(len(t2)):
    print(f"'{t2[i]}': {i+1},")
