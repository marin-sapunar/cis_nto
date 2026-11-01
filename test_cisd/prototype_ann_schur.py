import torch
import torch.nn as nn
import numpy as np
import time

# 1. Generate training data (realistic 4x4 W matrices)
np.random.seed(42)
n_samples = 50000
W = np.zeros((n_samples, 16))
targets = np.zeros(n_samples)
for i in range(n_samples):
    w = np.random.randn(4, 4) * 0.02
    w = (w + w.T) / 2
    W[i] = w.flatten()
    targets[i] = np.linalg.det(np.eye(4) + w)

# Train/test split
split = int(0.9 * n_samples)
X_train, X_test = torch.tensor(W[:split], dtype=torch.float32), torch.tensor(W[split:], dtype=torch.float32)
y_train, y_test = torch.tensor(targets[:split], dtype=torch.float32).unsqueeze(1), torch.tensor(targets[split:], dtype=torch.float32).unsqueeze(1)

# 2. Define the ANN
class SchurNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(16, 64),
            nn.SiLU(),
            nn.Linear(64, 64),
            nn.SiLU(),
            nn.Linear(64, 1)
        )
    def forward(self, x):
        return self.net(x)

model = SchurNet()
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
loss_fn = nn.MSELoss()

# 3. Training loop
print("Training ANN...")
t0 = time.time()
for epoch in range(500):
    model.train()
    optimizer.zero_grad()
    pred = model(X_train)
    loss = loss_fn(pred, y_train)
    loss.backward()
    optimizer.step()
    
    if (epoch + 1) % 100 == 0:
        print(f"Epoch {epoch+1}/500, MSE: {loss.item():.4e}")
t_train = time.time() - t0
print(f"Training time: {t_train:.2f} s")

# 4. Evaluate
model.eval()
with torch.no_grad():
    pred = model(X_test)
    abs_err = torch.abs(pred - y_test).numpy()
    rel_err = abs_err / np.abs(y_test.numpy())
    
print("\n=== Results ===")
print(f"Mean absolute error: {abs_err.mean():.4e}")
print(f"Max absolute error:  {abs_err.max():.4e}")
print(f"Mean relative error: {rel_err.mean():.4e}")
print(f"Max relative error:  {rel_err.max():.4e}")

# 5. Speed test
x_single = X_test[:1]
t0 = time.time()
with torch.no_grad():
    for _ in range(10000):
        _ = model(x_single)
t_ann = (time.time() - t0) / 10000 * 1e6

print(f"\nANN prediction time: {t_ann:.3f} µs per call")
