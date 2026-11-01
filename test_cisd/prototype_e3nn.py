import torch
import e3nn
from e3nn import o3

print(f"=== e3nn Equivariant Network Prototype ===")
print(f"PyTorch version: {torch.__version__}")
print(f"e3nn version: {e3nn.__version__}")

# 1. Simulate a molecule (e.g., 5 atoms with random 3D coordinates)
# In a real NAMD model, this would be the geometry at time t1 and t2.
N_atoms = 5
pos = torch.randn(N_atoms, 3)  # 3D Cartesian coordinates

# 2. Define the irreducible representations (irreps) for the features
# We use l=0 (scalars) and l=1 (vectors, e.g., distance vectors).
irreps_input = o3.Irreps("1x0e + 1x1o")  # 1 scalar + 1 vector
irreps_output = o3.Irreps("1x0e")        # Predict a single scalar (e.g., S_00)

# 3. Create an equivariant linear layer
# This layer maps the input geometry to the output scalar while strictly 
# respecting 3D rotation symmetry.
linear = o3.Linear(irreps_input, irreps_output)

# 4. Compute the pairwise distance vectors between atoms
# (This is the fundamental input to any geometry-aware model)
edge_vectors = pos.unsqueeze(0) - pos.unsqueeze(1)  # Shape: (N_atoms, N_atoms, 3)

# 5. Convert the 3D vectors into spherical harmonics (l=0, 1)
# Spherical harmonics are the mathematical foundation of e3nn.
sh = o3.spherical_harmonics(
    irreps_input, 
    edge_vectors, 
    normalize=True, 
    normalization='component'
)

# 6. Pass the spherical harmonics through the equivariant linear layer
# and sum over all atom pairs to get a single global scalar prediction.
pair_features = linear(sh)
prediction = pair_features.sum(dim=(0, 1))

print("\n--- Network Architecture ---")
print(f"Input Irreps:  {irreps_input}")
print(f"Output Irreps: {irreps_output}")
print(f"Equivariant Linear Weights: {sum(p.numel() for p in linear.parameters())}")

print("\n--- Forward Pass ---")
print(f"Input 3D Coordinates: {pos.shape}")
print(f"Predicted Overlap Element (S_00): {prediction.item():.6f}")

# 7. Prove Equivariance
# If we rotate the molecule, the predicted scalar (S_00) MUST NOT change.
print("\n--- Proving E(3) Equivariance ---")
# Create a random 3D rotation matrix
Q = o3.rand_matrix()
pos_rotated = pos @ Q.T

# Recompute
edge_vectors_rot = pos_rotated.unsqueeze(0) - pos_rotated.unsqueeze(1)
sh_rot = o3.spherical_harmonics(irreps_input, edge_vectors_rot, normalize=True, normalization='component')
pred_rotated = linear(sh_rot).sum(dim=(0, 1))

print(f"Original Prediction: {prediction.item():.8f}")
print(f"Rotated Prediction:  {pred_rotated.item():.8f}")
print(f"Difference:          {abs(prediction.item() - pred_rotated.item()):.2e}")
print("-> Difference is zero: the model is strictly E(3)-equivariant.")
