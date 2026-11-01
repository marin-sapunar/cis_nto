import numpy as np
import torch
import torch.nn as nn
import time

print("=== Physics-Constrained ML: Residual Learning for det(I+W) ===\n")

# ---- 1. Load real amyloid geometry ----
with open("real_data_amyloid.txt", "r") as f:
    header = f.readline().split()
    n_mo = int(header[0])
    S_mo = np.array([[float(x) for x in f.readline().split()] for _ in range(n_mo)])

# ---- 2. Build realistic W blocks ----
np.random.seed(42)
n_samples = 30000
W_data = np.zeros((n_samples, 16))
targets_det = np.zeros(n_samples)

n_occ = 20
A_ref = S_mo[:n_occ, :n_occ]
A_inv = np.linalg.inv(A_ref)

for i in range(n_samples):
    exc_r1, exc_r2 = np.random.choice(n_occ, 2, replace=False)
    exc_c1, exc_c2 = np.random.choice(n_occ, 2, replace=False)
    v1 = S_mo[exc_r1, :n_occ] - A_ref[exc_r1, :]
    v2 = S_mo[exc_r2, :n_occ] - A_ref[exc_r2, :]
    u3 = np.zeros(n_occ); u3[exc_c1] = 1.0
    u4 = np.zeros(n_occ); u4[exc_c2] = 1.0
    Av1, Av2 = A_inv @ v1, A_inv @ v2
    Au3, Au4 = A_inv @ u3, A_inv @ u4
    W = np.zeros((4, 4))
    W[0,0] = Av1[exc_r1]; W[0,1] = Av1[exc_r2]
    W[1,0] = Av2[exc_r1]; W[1,1] = Av2[exc_r2]
    W[0,2] = np.dot(v1, Au3); W[0,3] = np.dot(v1, Au4)
    W[1,2] = np.dot(v2, Au3); W[1,3] = np.dot(v2, Au4)
    W[2,0] = A_inv[exc_c1, exc_r1]; W[2,1] = A_inv[exc_c1, exc_r2]
    W[3,0] = A_inv[exc_c2, exc_r1]; W[3,1] = A_inv[exc_c2, exc_r2]
    W[2,2] = Au3[exc_c1]; W[2,3] = Au4[exc_c1]
    W[3,2] = Au3[exc_c2]; W[3,3] = Au4[exc_c2]
    W_data[i] = W.flatten()
    targets_det[i] = np.linalg.det(np.eye(4) + W)

# ---- 3. Feature augmentation ----
# Add the trace, Frobenius norm, and sum of squares as extra features
X = np.zeros((n_samples, 16 + 5))
X[:, :16] = W_data
X[:, 16] = np.trace(W_data.reshape(-1, 4, 4), axis1=1, axis2=2)
X[:, 17] = np.sum(W_data**2, axis=1)
X[:, 18] = np.linalg.norm(W_data.reshape(-1, 4, 4), axis=(1, 2))
X[:, 19] = targets_det - 1.0  # sanity feature (removed later)
X = X[:, :19]

# ---- 4. Train/test split ----
split = int(0.9 * n_samples)
X_train = torch.tensor(X[:split], dtype=torch.float32)
X_test = torch.tensor(X[split:], dtype=torch.float32)
y_det_train = torch.tensor(targets_det[:split], dtype=torch.float32).unsqueeze(1)
y_det_test = torch.tensor(targets_det[split:], dtype=torch.float32).unsqueeze(1)

# ---- 5. Define network ----
class SchurNet(nn.Module):
    def __init__(self, in_dim=19):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, 128), nn.SiLU(),
            nn.Linear(128, 128), nn.SiLU(),
            nn.Linear(128, 64), nn.SiLU(),
            nn.Linear(64, 1)
        )
    def forward(self, x): return self.net(x)

def train_and_eval(target_name, y_train, y_test, add_one=False):
    torch.manual_seed(42)
    model = SchurNet()
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = nn.MSELoss()

    t0 = time.time()
    for epoch in range(800):
        opt.zero_grad()
        pred = model(X_train)
        loss = loss_fn(pred, y_train)
        loss.backward()
        opt.step()
    t_train = time.time() - t0

    model.eval()
    with torch.no_grad():
        pred = model(X_test)
        if add_one:
            pred = pred + 1.0
        err = torch.abs(pred - y_test).numpy()
        rel = err / np.abs(y_test.numpy())
    return err.mean(), err.max(), rel.mean(), rel.max(), t_train

# ---- 6. Train three variants ----
print("Training Variant A: predict det(I+W) directly...")
errA = train_and_eval("A", y_det_train, y_det_test, add_one=False)
print("Training Variant B: predict delta = det(I+W) - 1 (residual)...")
delta_train = y_det_train - 1.0
delta_test = y_det_test - 1.0
errB = train_and_eval("B", delta_train, delta_test, add_one=True)
print("Training Variant C: predict log(det(I+W))...")
log_train = torch.log(y_det_train)
log_test = torch.log(y_det_test)
errC = train_and_eval("C", log_train, log_test, add_one=False)
# exp back
with torch.no_grad():
    model = SchurNet(); opt = torch.optim.Adam(model.parameters(), lr=1e-3); loss_fn = nn.MSELoss()
    for _ in range(800):
        opt.zero_grad(); loss = loss_fn(model(X_train), log_train); loss.backward(); opt.step()
    model.eval()
    predC = torch.exp(model(X_test))
    errC = (torch.abs(predC - y_det_test).numpy().mean(),
            torch.abs(predC - y_det_test).numpy().max(),
            (torch.abs(predC - y_det_test).numpy() / np.abs(y_det_test.numpy())).mean(),
            (torch.abs(predC - y_det_test).numpy() / np.abs(y_det_test.numpy())).max(),
            0.0)

# ---- 7. Report ----
print()
print(f"{'Variant':>10} | {'Mean abs':>12} | {'Max abs':>12} | {'Mean rel':>12} | {'Max rel':>12} | {'Train(s)':>8}")
print("-" * 82)
for name, r in [("Direct", errA), ("Residual", errB), ("Log", errC)]:
    print(f"{name:>10} | {r[0]:>12.4e} | {r[1]:>12.4e} | {r[2]:>12.4e} | {r[3]:>12.4e} | {r[4]:>8.2f}")

best = min([("Direct", errA), ("Residual", errB), ("Log", errC)], key=lambda x: x[1][0])
print(f"\n-> Best target: {best[0]}  (mean abs error {best[1][0]:.4e})")
print("-> Physics-constrained targets often beat direct regression because the")
print("   network only needs to learn the *correction* to the trivial det(I)=1.")
