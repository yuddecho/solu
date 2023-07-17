import matplotlib.pyplot as plt
import numpy as np

train_loss, train_acc = [], []
chang_loss, chang_acc = [], []
nesg_loss, nesg_acc = [], []

data_file = '../../data/loos.csv'
with open(data_file, 'r', encoding='utf-8') as r:
    lines = r.readlines()
    for line in lines:
        line = line.strip()

        if not line:
            break

        # 分割数据
        line_list = line.split(',')
        temp = []
        for item in line_list:
            temp.append(round(float(item), 3))

        line_list = temp

        train_loss.append(line_list[0])
        chang_loss.append(line_list[1])
        nesg_loss.append(line_list[2])

        train_acc.append(line_list[3])
        chang_acc.append(line_list[4])
        nesg_acc.append(line_list[5])

# 规整
y_max, y_min = max(chang_loss), min(chang_loss)
temp = []
for item in chang_loss:
    temp.append(item / (y_max - y_min))
chang_loss = temp

y_max, y_min = max(nesg_loss), min(nesg_loss)
temp = []
for item in nesg_loss:
    temp.append(item / (y_max - y_min))
nesg_loss = temp

plt.figure(figsize=(14, 6))

x = range(len(train_loss))
print(train_loss)

# 图一
plt.subplot(121)

plt.scatter(x, train_loss, s=2, label='tt-pdb')
plt.scatter(x, chang_loss, s=4, label='chang')
plt.scatter(x, nesg_loss, s=2, label='nesg')

plt.legend()

new_ticks = np.linspace(0, 1, 5)
plt.yticks(new_ticks)

plt.title('Change of loss function value per epoch')
plt.xlabel('epoch')
plt.ylabel('loss')


# 图二
plt.subplot(122)

plt.scatter(x, train_acc, s=2, label='tt-pdb')
plt.scatter(x, chang_acc, s=4, label='chang')
plt.scatter(x, nesg_acc, s=2, label='nesg')

plt.legend()

new_ticks = np.linspace(0, 1, 5)
plt.yticks(new_ticks)

plt.title('Change of acc per epoch')
plt.xlabel('epoch')
plt.ylabel('acc')

plt.savefig('../../data/res.png')

plt.show()
