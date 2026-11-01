import numpy as np
import torch
import torch.nn as nn
import time

print("=== FINAL ML PROTOTYPE: Real Amyloid Schur Blocks ===")

# 1. Load the REAL amyloid MO overlap matrix
with open("real_data_amyloid.txt", "r") as f:
    header = f.readline().split()
    n_mo = int(header[0])
    S_mo = np.zeros((n_mo, n_mo))
    for i in range(n_mo):
        S_mo[i, :] = np.array([float(x) for x in f.readline().split()])
print(f"Loaded real MO overlap matrix: {n_mo}x{n_mo}")

# 2. Extract realistic W blocks
np.random.seed(42)
n_samples = 30000
W_data = np.zeros((n_samples, 16))
targets = np.zeros(n_samples)

n_occ = 20
A_ref = S_mo[:n_occ, :n_occ]
A_inv = np.linalg.inv(A_ref)

for i in range(n_samples):
    exc_r1, exc_r2 = np.random.choice(n_occ, 2, replace=False)
    exc_c1, exc_c2 = np.random.choice(n_occ, 2, replace=False)
    
    v1 = S_mo[exc_r1, :n_occ] - A_ref[exc_r1, :]
    v2 = S_mo[exc_r2, :n_occ] - A_ref[exc_r2, :]
    
    u3 = np.zeros(n_occ)
    u4 = np.zeros(n_occ)
    u3[exc_c1] = 1.0
    u4[exc_c2] = 1.0
    
    Ainv_v1 = A_inv @ v1
    Ainv_v2 = A_inv @ v2
    Ainv_u3 = A_inv @ u3
    Ainv_u4 = A_inv @ u4
    
    W = np.zeros((4, 4))
    W[0, 0] = Ainv_v1[exc_r1]
    W[0, 1] = Ainv_v1[exc_r2]
    W[1, 0] = Ainv_v2[exc_r1]
    W[1, 1] = Ainv_v2[exc_r2]
    
    W[0, 2] = np.dot(v1, Ainv_u3)
    W[0, 3] = np.dot(v1, Ainv_u4)
    W[1, 2] = np.dot(v2, Ainv_u3)
    W[1, 3] = np.dot(v2, Ainv_u4)
    
    W[2, 0] = A_inv[exc_c1, exc_r1]
    W[2, 1] = A_inv[exc_c1, exc_r2]
    W[3, 0] = A_inv[exc_c2, exc_r1]
    W[3, 1] = A_inv[exc_c2, exc_r2]
    
    W[2, 2] = Ainv_u3[exc_c1]
    W[2, 3] = Ainv_u4[exc_c1]
    W[3, 2] = Ainv_u3[exc_c2]
    W[3, 3] = Ainv_u4[exc_c2]
    
    W_data[i, :] = W.flatten()
    targets[i] = np.linalg.det(np.eye(4) + W)

# 3. Train/Test Split
split = int(0.9 * n_samples)
X_train = torch.tensor(W_data[:split], dtype=torch.float32)
y_train = torch.tensor(targets[:split], dtype=torch.float32).unsqueeze(1)
X_test = torch.tensor(W_data[split:], dtype=torch.float32)
y_test = torch.tensor(targets[split:], dtype=torch.float32).unsqueeze(1)

# 4. Define the Final ANN
class FinalSchurNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(16, 128), nn.SiLU(),
            nn.Linear(128, 128), nn.SiLU(),
            nn.Linear(128, 64), nn.SiLU(),
            nn.Linear(64, 1)
        )
    def forward(self, x): return self.net(x)

model = FinalSchurNet()
opt = torch.optim.Adam(model.parameters(), lr=1e-3)
loss_fn = nn.MSELoss()

print("Training ANN on real amyloid data...")
t0 = time.time()
for epoch in range(500):
    opt.zero_grad()
    loss = loss_fn(model(X_train), y_train)
    loss.backward()
    opt.step()
t_train = time.time() - t0

model.eval()
with torch.no_grad():
    pred = model(X_test)
    abs_err = torch.abs(pred - y_test).numpy()
    rel_err = abs_err / np.abs(y_test.numpy())

print(f"\n=== Final Results ===")
print(f"Train time: {t_train:.2f} s")
print(f"Mean absolute error: {abs_err.mean():.4e}")
print(f"Max absolute error:  {abs_err.max():.4e}")
print(f"Mean relative error: {rel_err.mean():.4e}")
print(f"Max relative error:  {rel_err.max():.4e}")
