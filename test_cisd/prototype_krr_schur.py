import numpy as np
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import train_test_split
import time

np.random.seed(42)
n_samples = 50000

# Generate realistic 4x4 W matrices
W = np.zeros((n_samples, 16))
targets = np.zeros(n_samples)
for i in range(n_samples):
    w = np.random.randn(4, 4) * 0.02
    w = (w + w.T) / 2
    W[i] = w.flatten()
    targets[i] = np.linalg.det(np.eye(4) + w)

# Split off the test set FIRST (from all 50k)
W_train_full, X_test, y_train_full, y_test = train_test_split(
    W, targets, test_size=0.1, random_state=42)

# Then subsample the training set to 2,000 (KRR memory limit)
n_train_sub = 2000
idx = np.random.choice(len(W_train_full), n_train_sub, replace=False)
X_train = W_train_full[idx]
y_train = y_train_full[idx]

print(f"KRR training set: {n_train_sub} samples (subsampled from 45,000)")
print(f"KRR test set:     {len(X_test)} samples")
print("Training KRR...")

t0 = time.time()
krr = KernelRidge(kernel='rbf', alpha=1e-8, gamma=1.0)
krr.fit(X_train, y_train)
t_train = time.time() - t0

y_pred = krr.predict(X_test)
abs_err = np.abs(y_pred - y_test)
rel_err = abs_err / np.abs(y_test)

print(f"\n=== Results ===")
print(f"Train time: {t_train:.3f} s")
print(f"Mean absolute error: {abs_err.mean():.4e}")
print(f"Max absolute error:  {abs_err.max():.4e}")
print(f"Mean relative error: {rel_err.mean():.4e}")
print(f"Max relative error:  {rel_err.max():.4e}")
print(f"R^2 score:            {krr.score(X_test, y_test):.6f}")

# Speed test (KRR prediction scales with training size)
t0 = time.time()
for i in range(1000):
    _ = krr.predict(X_test[:1])
t_krr = (time.time() - t0) / 1000 * 1e6
print(f"\nKRR inference time: {t_krr:.2f} µs per call")
print(f"Note: KRR inference cost scales linearly with training set size.")
print(f"      ANN inference was 23.5 µs (constant).")
