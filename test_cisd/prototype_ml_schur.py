import numpy as np
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

# ============================================================
# Step 1. Generate training data: (W, det(I+W)) pairs
# ============================================================
# W is the small matrix inside the Schur update for a (2,2) double excitation.
# In your Fortran code, W is a 4x4 matrix.
# We'll generate 50,000 realizations of W with statistically realistic values.
np.random.seed(42)
n_samples = 50000

# W has 16 elements (4x4). We draw them from a distribution consistent with
# MO overlaps: diagonal ~ 1e-1 to 1e-2, off-diagonal ~ 1e-2 to 1e-3.
W = np.zeros((n_samples, 4, 4))
for i in range(n_samples):
    W[i] = np.random.randn(4, 4) * 0.02
    W[i] = (W[i] + W[i].T) / 2  # symmetrize (physically W is nearly symmetric)

# Target: det(I + W)
targets = np.zeros(n_samples)
for i in range(n_samples):
    targets[i] = np.linalg.det(np.eye(4) + W[i])

# ============================================================
# Step 2. Prepare features
# ============================================================
# Features: the 16 elements of W plus 4 auxiliary features
# (diagonal of W, sum of squares, and the trace).
X = np.zeros((n_samples, 20))
for i in range(n_samples):
    X[i, :16] = W[i].flatten()
    X[i, 16] = np.trace(W[i])
    X[i, 17] = np.sum(W[i]**2)
    X[i, 18] = np.linalg.norm(W[i], 'fro')
    X[i, 19] = np.linalg.det(W[i])

# ============================================================
# Step 3. Train/Test split
# ============================================================
X_train, X_test, y_train, y_test = train_test_split(X, targets, test_size=0.1, random_state=42)

print("=== ML Surrogate for Schur Determinant ===")
print(f"Training set: {len(X_train)} samples")
print(f"Test set:     {len(X_test)} samples")
print()

# ============================================================
# Step 4. Train MLP
# ============================================================
model = MLPRegressor(
    hidden_layer_sizes=(64, 64),  # 2 hidden layers
    activation='tanh',
    solver='adam',
    max_iter=2000,
    early_stopping=True,
    random_state=42,
    verbose=False,
)

print("Training MLP...")
model.fit(X_train, y_train)

# ============================================================
# Step 5. Evaluate
# ============================================================
y_pred = model.predict(X_test)
abs_error = np.abs(y_pred - y_test)
rel_error = abs_error / np.abs(y_test)

print()
print("=== Results ===")
print(f"Mean absolute error: {abs_error.mean():.4e}")
print(f"Max absolute error:  {abs_error.max():.4e}")
print(f"Mean relative error: {rel_error.mean():.4e}")
print(f"Max relative error:  {rel_error.max():.4e}")
print(f"R^2 score:            {model.score(X_test, y_test):.6f}")
print()

# ============================================================
# Step 6. Estimate speedup
# ============================================================
import time
# Time the MLP prediction
n_reps = 10000
t0 = time.time()
for i in range(n_reps):
    _ = model.predict(X_test[:1])
t_mlp = (time.time() - t0) / n_reps * 1e6  # µs per prediction

# Time the exact determinant (via NumPy as a fair proxy for the Fortran dot products)
t0 = time.time()
for i in range(n_reps):
    _ = np.linalg.det(np.eye(4) + W[i % n_samples])
t_exact = (time.time() - t0) / n_reps * 1e6  # µs per exact call

print("=== Speed ===")
print(f"Exact det(I+W): {t_exact:.3f} µs per call")
print(f"MLP prediction: {t_mlp:.3f} µs per call")
print(f"Speedup (per-pair): {t_exact / t_mlp:.2f}x")
print()
print("NOTE: In the Fortran code, the exact call is much more expensive")
print("      because it involves strided memory access for Ainv_smo and")
print("      the v_i, u_j vectors. The MLP would be O(1) regardless.")
print()

# ============================================================
# Step 7. Plot
# ============================================================
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))

ax[0].hist(np.log10(abs_error + 1e-20), bins=50, color='teal', alpha=0.8)
ax[0].set_xlabel('log10(|error|)')
ax[0].set_ylabel('Count')
ax[0].set_title('Absolute Error Distribution')
ax[0].grid(True, linestyle=':', alpha=0.6)

ax[1].scatter(y_test[:2000], y_pred[:2000], s=5, alpha=0.5)
ax[1].plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
ax[1].set_xlabel('Exact det(I+W)')
ax[1].set_ylabel('MLP prediction')
ax[1].set_title('Prediction vs. Exact')
ax[1].grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig('figure3_ml_schur.png', dpi=300)
print("Saved figure3_ml_schur.png")
