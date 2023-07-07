from sklearn.datasets import make_multilabel_classification
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from sklearn.model_selection import train_test_split

X, y = make_multilabel_classification(n_samples=1000,
                                      n_features=10,
                                      n_classes=2,
                                      n_labels=2,
                                      random_state=1)

n_inputs, n_outputs = X.shape[1], y.shape[1]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)
X_train = torch.from_numpy(X_train).float()
X_test = torch.from_numpy(X_test).float()
y_train = torch.from_numpy(y_train).float()
y_test = torch.from_numpy(y_test).float()

train_data = [(X, y) for X, y in zip(X_train, y_train)]
train_loader = DataLoader(train_data, batch_size=64, shuffle=True)


class MLP(nn.Module):
    def __init__(self, n_inputs, n_outputs, num_hiddens):
        super(MLP, self).__init__()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(n_inputs, num_hiddens),
            nn.ReLU(),
            nn.Linear(num_hiddens, n_outputs),
            nn.Sigmoid())

    def forward(self, x):
        outputs = self.linear_relu_stack(x)
        return outputs


num_hiddens = 30
model = MLP(n_inputs, n_outputs, num_hiddens)
print(model)

loss = nn.BCELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)


def evaluate_accuracy(X, y, model):
    pred = model(X)
    correct = sum(row.all().int().item() for row in (pred.ge(0.5) == y))
    n = y.shape[0]
    return correct / n


def train(train_loader, X_test, y_test, model, loss, num_epochs, batch_size,
          optimizer):
    batch_count = 0
    for epoch in range(num_epochs):
        train_l_sum, train_acc_sum, n = 0.0, 0.0, 0
        for X, y in train_loader:
            pred = model(X)
            l = loss(pred, y)
            optimizer.zero_grad()
            l.backward()
            optimizer.step()
            train_l_sum += l.item()
            train_acc_sum += sum(row.all().int().item()
                                 for row in (pred.ge(0.5) == y))
            n += y.shape[0]
            batch_count += 1
        test_acc = evaluate_accuracy(X_test, y_test, model)
        print(
            'epoch %d, loss %.4f, train acc %.3f, test acc %.3f'
            % (epoch + 1, train_l_sum / batch_count, train_acc_sum / n,
               test_acc))


num_epochs, batch_size = 20, 64
train(train_loader, X_test, y_test, model, loss, num_epochs, batch_size, optimizer)
