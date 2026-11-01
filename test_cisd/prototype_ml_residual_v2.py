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

print(f"Loaded {n_samples} samples.")
print(f"  det(I+W) range: [{targets_det.min():.4f}, {targets_det.max():.4f}]")
print(f"  Fraction with det < 0: {(targets_det < 0).mean()*100:.1f}%")
print(f"  Fraction with det ≈ 1 (within 1%): {(np.abs(targets_det - 1) < 0.01).mean()*100:.1f}%\n")

# ---- 3. Feature construction ----
W4 = W_data.reshape(-1, 4, 4)
X = np.zeros((n_samples, 19))
X[:, :16] = W_data
X[:, 16] = np.trace(W4, axis1=1, axis2=2)
X[:, 17] = np.sum(W4**2, axis=(1, 2))
X[:, 18] = np.linalg.norm(W4, axis=(1, 2))

split = int(0.9 * n_samples)
X_train = torch.tensor(X[:split], dtype=torch.float32)
X_test = torch.tensor(X[split:], dtype=torch.float32)
y_det_test = torch.tensor(targets_det[split:], dtype=torch.float32).unsqueeze(1)

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

def train_variant(name, target_train, target_test, transform_fn, inverse_fn):
    torch.manual_seed(42)
    model = SchurNet()
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = nn.MSELoss()
    
    # Apply transformation to targets
    y_train_t = transform_fn(target_train)
    y_test_t = transform_fn(target_test)
    
    t0 = time.time()
    for epoch in range(800):
        opt.zero_grad()
        loss = loss_fn(model(X_train), y_train_t)
        loss.backward()
        opt.step()
    t_train = time.time() - t0
    
    model.eval()
    with torch.no_grad():
        pred_t = model(X_test)
        pred = inverse_fn(pred_t)
        err = torch.abs(pred - y_det_test).numpy()
        rel = err / np.abs(y_det_test.numpy())
    return err.mean(), err.max(), rel.mean(), rel.max(), t_train, pred.numpy().flatten(), y_det_test.numpy().flatten()

y_train_det = torch.tensor(targets_det[:split], dtype=torch.float32).unsqueeze(1)
y_test_det = torch.tensor(targets_det[split:], dtype=torch.float32).unsqueeze(1)

# Variant A: direct
resA = train_variant("A", y_train_det, y_test_det,
                     transform_fn=lambda x: x,
                     inverse_fn=lambda x: x)

# Variant B: residual
resB = train_variant("B", y_train_det, y_test_det,
                     transform_fn=lambda x: x - 1.0,
                     inverse_fn=lambda x: x + 1.0)

# Variant C: log with sign handling (signed_log)
def signed_log(x):
    return torch.sign(x) * torch.log(torch.abs(x) + 1e-12)
def inverse_signed_log(x):
    return torch.sign(x) * (torch.exp(torch.abs(x)) - 1e-12)
resC = train_variant("C", y_train_det, y_test_det,
                     transform_fn=signed_log,
                     inverse_fn=inverse_signed_log)

# Report
print(f"{'Variant':>10} | {'Mean abs':>12} | {'Max abs':>12} | {'Mean rel':>12} | {'Max rel':>12} | {'Train(s)':>8}")
print("-" * 85)
for name, r in [("Direct", resA), ("Residual", resB), ("Signed log", resC)]:
    print(f"{name:>10} | {r[0]:>12.4e} | {r[1]:>12.4e} | {r[2]:>12.4e} | {r[3]:>12.4e} | {r[4]:>8.2f}")

best = min([("Direct", resA), ("Residual", resB), ("Signed log", resC)], key=lambda x: x[1][0])
print(f"\n-> Best target: {best[0]}  (mean abs error {best[1][0]:.4e})")
