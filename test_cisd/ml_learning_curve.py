import numpy as np
import time
import torch
import torch.nn as nn
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

np.random.seed(42)
n_samples = 50000
W = np.zeros((n_samples, 16))
targets = np.zeros(n_samples)
for i in range(n_samples):
    w = np.random.randn(4, 4) * 0.02
    w = (w + w.T) / 2
    W[i] = w.flatten()
    targets[i] = np.linalg.det(np.eye(4) + w)

X_pool, X_test, y_pool, y_test = train_test_split(W, targets, test_size=10000, random_state=42)
X_test_t = torch.tensor(X_test, dtype=torch.float32)
y_test_t = torch.tensor(y_test, dtype=torch.float32).unsqueeze(1)

class SchurNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(16, 64), nn.SiLU(),
            nn.Linear(64, 64), nn.SiLU(),
            nn.Linear(64, 1)
        )
    def forward(self, x): return self.net(x)

train_sizes = [100, 500, 1000, 2000, 5000, 10000, 20000, 40000]
krr_errors, ann_errors = [], []

print(f"{'Size':>6} | {'KRR Error':>12} | {'ANN Error':>12} | {'KRR Time':>10} | {'ANN Time':>10}")
print("-" * 65)

FIXED_EPOCHS = 800  # same for every training size

for size in train_sizes:
    idx = np.random.choice(len(X_pool), size, replace=False)
    X_tr, y_tr = X_pool[idx], y_pool[idx]
    X_tr_t = torch.tensor(X_tr, dtype=torch.float32)
    y_tr_t = torch.tensor(y_tr, dtype=torch.float32).unsqueeze(1)

    # --- KRR ---
    t0 = time.time()
    krr = KernelRidge(kernel='rbf', alpha=1e-8, gamma=1.0)
    krr.fit(X_tr, y_tr)
    k_time = time.time() - t0
    k_err = np.abs(krr.predict(X_test) - y_test).mean()

    # --- ANN (fixed epochs) ---
    t0 = time.time()
    model = SchurNet()
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = nn.MSELoss()
    for _ in range(FIXED_EPOCHS):
        opt.zero_grad()
        loss = loss_fn(model(X_tr_t), y_tr_t)
        loss.backward()
        opt.step()
    a_time = time.time() - t0
    model.eval()
    with torch.no_grad():
        a_err = torch.abs(model(X_test_t) - y_test_t).mean().item()

    krr_errors.append(k_err); ann_errors.append(a_err)
    print(f"{size:>6} | {k_err:>12.2e} | {a_err:>12.2e} | {k_time:>9.2f}s | {a_time:>9.2f}s")

# Plot
fig, ax1 = plt.subplots(figsize=(8, 5))
ax1.loglog(train_sizes, krr_errors, 'o-', color='blue', label='KRR Error')
ax1.loglog(train_sizes, ann_errors, 's-', color='red', label='ANN Error')
ax1.set_xlabel('Training Set Size')
ax1.set_ylabel('Mean Absolute Error')
ax1.set_title('Learning Curve: KRR vs ANN for Schur Determinant')
ax1.grid(True, which="both", ls="--", alpha=0.5)
ax1.legend()
plt.tight_layout()
plt.savefig('figure4_ml_learning_curve.png', dpi=300)
print("\nSaved figure4_ml_learning_curve.png")
