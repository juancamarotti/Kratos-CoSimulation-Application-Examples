import KratosMultiphysics as KM
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
import numpy as np 
import math

# ------------------- Load Model Parts -------------------

model = KM.Model()
source_model_part = model.CreateModelPart("Source")
target_model_part = model.CreateModelPart("Target")

KM.ModelPartIO("plate_fine").ReadModelPart(source_model_part)
KM.ModelPartIO("plate_coarse").ReadModelPart(target_model_part)

# ------------------- RBF Mapping Function -------------------

def rbf_map_field(source_model_part, target_model_part, analytical_function, rbf_function, variable=KM.TEMPERATURE, support_radius=None):
    """
    Maps an analytical field from a source model part to a target using user-defined RBF interpolation.

    Parameters:
    - source_model_part: Kratos ModelPart (fine mesh)
    - target_model_part: Kratos ModelPart (coarse mesh)
    - analytical_function: function f(x, y, z)
    - rbf_function: function phi(r, R) -> float (compact support should use R)
    - variable: Kratos variable (e.g., TEMPERATURE)
    - support_radius: Optional fixed support radius (for compact RBFs)
    """

    source_coords = []
    source_values = []

    for node in source_model_part.Nodes:
        x, y, z = node.X, node.Y, node.Z
        source_coords.append([x, y, z])
        value = analytical_function(x, y, z)
        source_values.append(value)
        node.SetValue(variable, value)

    source_coords = np.array(source_coords)
    source_values = np.array(source_values)

    # --- Estimate support radius if needed ---
    def estimate_support_radius(coords, k=2.5):
        coords = np.array(coords)
        total_dist = 0.0
        count = 0
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                total_dist += np.linalg.norm(coords[i] - coords[j])
                count += 1
        avg_spacing = total_dist / count
        return k * avg_spacing

    R = support_radius or estimate_support_radius(source_coords)

    # --- Build interpolation matrix ---
    def build_rbf_matrix(coords):
        n = len(coords)
        U = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                r = np.linalg.norm(coords[i] - coords[j])
                U[i, j] = rbf_function(r, R)
        return U

    def solve_coeffs(coords, values):
        A = build_rbf_matrix(coords)
        print(f"[INFO] RBF Matrix rank: {np.linalg.matrix_rank(A)}")
        print(f"[INFO] RBF Matrix condition number: {np.linalg.cond(A)}")
        try:
            x = np.linalg.solve(A, values)
        except np.linalg.LinAlgError as e:
            raise RuntimeError(f"[ERROR] Singular matrix: {e}")
        return x

    c = solve_coeffs(source_coords, source_values)

    # --- Interpolate at target nodes ---
    def evaluate_rbf_at(x, coords, c):
        rbf_sum = 0.0
        for cj, xj in zip(c, coords):
            r = np.linalg.norm(np.array(x) - xj)
            rbf_sum += cj * rbf_function(r, R)
        return rbf_sum

    for node in target_model_part.Nodes:
        x, y, z = node.X, node.Y, node.Z
        value = evaluate_rbf_at([x, y, z], source_coords, c)
        node.SetValue(variable, value)

# ------------------- Analytical Field -------------------

field = lambda x, y, z: math.sin(math.pi * (x))  # Shifted to avoid symmetry

# ------------------- Compact RBF Function -------------------

def wendland_c2_rbf(r, R):
    q = r / R
    if q >= 1.0:
        return 0.0
    return (1 - q)**4 * (4 * q + 1)

# ------------------- Multiquadric -------------------

def multiquadric_rbf(r, R=1e-3):
    return np.sqrt(r**2 + R**2)

# ------------------- Thin-plate Spline -------------------

def tps_rbf(r, R=1e-3):
    return r**2 * np.log(r**2) if r > 0 else 0.0

# ------------------- Gaussians -------------------

def gaussian_rbf(r, R=1e-3):
    return np.exp(-(R*r)**2)

# ------------------- Run Mapping -------------------

rbf_map_field(
    source_model_part,
    target_model_part,
    analytical_function=field,
    rbf_function=gaussian_rbf,
    variable=KM.TEMPERATURE
)

# ------------------- VTK Output -------------------

vtk_parameters = KM.Parameters(r"""
{
    "model_part_name"              : "Target",
    "output_control_type"          : "step",
    "output_interval"              : 1,
    "output_precision"             : 6,
    "output_sub_model_parts"       : false,
    "output_path"                  : "vtk_output",
    "save_output_files_in_folder"  : true,
    "nodal_data_value_variables" : ["TEMPERATURE"]
}
""")

vtk_output = VtkOutputProcess(model, vtk_parameters)
vtk_output.PrintOutput()  # One-shot write
