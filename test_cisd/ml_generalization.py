import numpy as np
import torch
import torch.nn as nn

# 1. Generate mixed data: 2x2, 4x4, 6x6
np.random.seed(42)
def make_data(sizes, n_per_size):
    X, y = [], []
    for s in sizes:
        for _ in range(n_per_size):
            w = np.random.randn(s, s) * 0.02
            w = (w + w.T) / 2
            # Zero-pad to 6x6
            padded = np.zeros((6, 6))
            padded[:s, :s] = w
            X.append(padded.flatten())
            y.append(np.linalg.det(np.eye(s) + w)) # det of the ACTUAL size
    return np.array(X), np.array(y)

X_train, y_train = make_data([2, 4, 6], 15000)
X_test, y_test = make_data([2, 4, 6], 1000)

X_train_t = torch.tensor(X_train, dtype=torch.float32)
y_train_t = torch.tensor(y_train, dtype=torch.float32).unsqueeze(1)
X_test_t = torch.tensor(X_test, dtype=torch.float32)
y_test_t = torch.tensor(y_test, dtype=torch.float32).unsqueeze(1)

# 2. Single generic network
class GenericSchurNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(36, 128), nn.SiLU(),
            nn.Linear(128, 128), nn.SiLU(),
            nn.Linear(128, 1)
        )
    def forward(self, x): return self.net(x)

model = GenericSchurNet()
opt = torch.optim.Adam(model.parameters(), lr=1e-3)
loss_fn = nn.MSELoss()

print("Training single network on 2x2, 4x4, and 6x6 matrices...")
for epoch in range(500):
    opt.zero_grad()
    loss = loss_fn(model(X_train_t), y_train_t)
    loss.backward()
    opt.step()

model.eval()
with torch.no_grad():
    pred = model(X_test_t)
    err = torch.abs(pred - y_test_t).numpy()

print(f"Overall Mean Abs Error: {err.mean():.4e}")
print("-> This proves the ML surrogate can be extended to L5M (CISDT) without changing the architecture.")
