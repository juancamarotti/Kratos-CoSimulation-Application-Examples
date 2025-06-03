import KratosMultiphysics as KM
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
import numpy as np 
import math

# ------------------- Load Model Parts -------------------

model = KM.Model()
source_model_part = model.CreateModelPart("Source")
target_model_part = model.CreateModelPart("Target")

KM.ModelPartIO("plate_source").ReadModelPart(source_model_part)
KM.ModelPartIO("plate_target").ReadModelPart(target_model_part)

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
        d = 0  # 3D case: x, y, z
        U = np.zeros((n, n))
        P = np.ones((n, d + 1))  # [1, x, y, z]

        # Fill RBF matrix and polynomial part
        for i in range(n):
            #P[i, 1:] = coords[i]
            for j in range(n):
                r = np.linalg.norm(coords[i] - coords[j])
                U[i, j] = rbf_function(r, R)

        print(f"[INFO] P Matrix rank: {np.linalg.matrix_rank(P)}")

        # Augmented system
        top = np.hstack((U, P))
        bottom = np.hstack((P.T, np.zeros((d + 1, d + 1))))
        A = np.vstack((top, bottom))

        return A, n, d  # also return 'n' to split the solution vector later

    def solve_coeffs(coords, values):
        A, n, poly_dim= build_rbf_matrix(coords)
        print(f"[INFO] RBF Matrix rank: {np.linalg.matrix_rank(A)}")
        print(f"[INFO] RBF Matrix condition number: {np.linalg.cond(A)}")

        # ConstruCT extended system
        b = np.concatenate((values, np.zeros(poly_dim+1)))
        
        try:
            x = np.linalg.solve(A, b)
        except np.linalg.LinAlgError as e:
            raise RuntimeError(f"[ERROR] Singular matrix: {e}")
        return x[:n], x[n:] # c, poly_coeffs

    c, poly_coeffs = solve_coeffs(source_coords, source_values)

    # --- Interpolate at target nodes ---
    def evaluate_rbf_at(x, coords, c, poly_coeffs):
        """
        Evaluates the RBF interpolant at a point x = [x, y, z].

        Parameters:
        - x: list or array of coordinates [x, y, z]
        - coords: list of source node coordinates
        - c: RBF weights (from solving the system)
        - poly_coeffs: polynomial coefficients (can be just [b0] or [b0, b1, b2, b3])

        Returns:
        - The interpolated value at point x
        """

        # RBF sum from each node
        rbf_sum = 0.0
        for cj, xj in zip(c, coords):
            r = np.linalg.norm(np.array(x) - xj)
            rbf_sum += cj * rbf_function(r, R)

        # Polynomial part: [1] or [1, x, y, z]
        if len(poly_coeffs) == 1:
            poly = poly_coeffs[0]
        else:
            b0, b1, b2, b3 = poly_coeffs
            poly = b0 + b1 * x[0] + b2 * x[1] + b3 * x[2]

        return rbf_sum + poly

    for node in target_model_part.Nodes:
        x, y, z = node.X, node.Y, node.Z
        value = evaluate_rbf_at([x, y, z], source_coords, c, poly_coeffs)
        node.SetValue(variable, value)

    def compute_rmse(target_model_part, analytical_function, variable):
        errors = []
        for node in target_model_part.Nodes:
            x, y, z = node.X, node.Y, node.Z
            true_value = analytical_function(x, y, z)
            predicted_value = node.GetValue(variable)
            errors.append((true_value - predicted_value) ** 2)

        rmse = math.sqrt(sum(errors) / len(errors))
        print(f"[INFO] RMSE = {rmse:.6e}")
        return rmse

    rmse = compute_rmse(target_model_part, analytical_function, variable)

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
    rbf_function=tps_rbf,
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
