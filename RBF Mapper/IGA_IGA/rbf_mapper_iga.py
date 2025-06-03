# Kratos Imports 
import KratosMultiphysics as Kratos
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA

# External imports
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.tri as tri

# ------------------- RBF Functions -------------------
def wendland_c2_rbf(r, R):
    q = r / R
    if q >= 1.0:
        return 0.0
    return (1 - q) ** 4 * (4 * q + 1)

def multiquadric_rbf(r, R=1e-3):
    return np.sqrt(r**2 + R**2)

def tps_rbf(r, R=1e-3):
    return r**2 * np.log(r**2) if r > 0 else 0.0

def gaussian_rbf(r, R=1e-3):
    return np.exp(-(R * r) ** 2)

# ------------------- Analytical field ----------------------------
def analytical_field(x, y, z):
    return math.sin(math.pi * x) + x**3 # Non-symmetric field

def analytical_force(x, y, z):
    return x**2 + np.sin(x)

# ------------------- HELPER FUNCTIONS -------------------

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

def build_origin_rbf_matrix(coords, rbf_function, R):
    n = len(coords)
    d = 0  # 3D: x, y, z
    PHI = np.zeros((n, n))
    Q = np.ones((n, d + 1))

    for i in range(n):
        for j in range(n):
            r = np.linalg.norm(coords[i] - coords[j])
            PHI[i, j] = rbf_function(r, R)

    top = np.hstack((PHI, Q))
    bottom = np.hstack((Q.T, np.zeros((d + 1, d + 1))))
    A = np.vstack((top, bottom))
    return A, n, d

def build_origin_destination_rbf_matrix(origin_coords, destination_coords, rbf_function, R):
    n_dest = len(destination_coords)
    n_origin = len(origin_coords)
    d = 0  # 3D: x, y, z
    PHI = np.zeros((n_dest, n_origin))
    Q = np.ones((n_dest, d + 1))

    for i in range(n_dest):
        for j in range(n_origin):
            r = np.linalg.norm(destination_coords[i] - origin_coords[j])
            PHI[i, j] = rbf_function(r, R)

    C = np.hstack((PHI, Q))
    
    return C

def compute_rmse(destination_model_part, analytical_function, variable):
    errors = []
    for element in destination_model_part.Elements:
        element_geometry = element.GetGeometry()
        x, y, z = element_geometry.Center().X, element_geometry.Center().Y, element_geometry.Center().Z
        true_value = analytical_function(x, y, z)
        predicted_value = element.GetValue(variable)
        errors.append((true_value - predicted_value) ** 2)

    rmse = math.sqrt(sum(errors) / len(errors))
    print(f"[INFO] RMSE = {rmse:.6e}")
    return rmse

def compute_total_force(model_part, variable):
    total_force = np.zeros(3)
    for element in model_part.Elements:
        force = element.GetValue(variable)
        if hasattr(force, '__len__') and len(force) == 3:
            total_force += np.array(force)
        else:
            total_force += np.array([force, 0, 0])  # Asumimos dirección X si es escalar
    return total_force

def tripcolor_plot(model_part, variable, title="Field Plot", cmap='viridis', filename=None):
    coords = []
    values = []

    for element in model_part.Elements:
        geom = element.GetGeometry()
        x, y = geom.Center().X, geom.Center().Y
        val = element.GetValue(variable)
        
        # Handle scalar or vector (e.g. POINT_LOAD)
        if hasattr(val, '__len__') and len(val) == 3:
            val = np.linalg.norm(val)
        else:
            val = float(val)

        coords.append([x, y])
        values.append(val)

    coords = np.array(coords)
    values = np.array(values)

    # Create triangulation
    triangulation = tri.Triangulation(coords[:, 0], coords[:, 1])

    # Plot with tripcolor
    plt.figure(figsize=(8, 6))
    tpc = plt.tripcolor(triangulation, values, shading='gouraud', cmap=cmap)
    plt.colorbar(tpc, label=variable.Name())
    plt.title(title)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis("equal")
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300)
        plt.close()
    else:
        plt.show()

def tripcolor_error_plot(destination_model_part, analytical_function, variable, title="Error Field Plot", cmap='inferno', filename=None):
    coords = []
    errors = []
    
    for element in destination_model_part.Elements:
        geom = element.GetGeometry()
        x, y, z = geom.Center().X, geom.Center().Y, geom.Center().Z
        predicted = element.GetValue(variable)

        # Handle scalar or vector variable
        if hasattr(predicted, '__len__') and len(predicted) == 3:
            predicted = np.linalg.norm(predicted)
        else:
            predicted = float(predicted)

        true_value = analytical_function(x, y, z)
        error = abs(predicted - true_value)

        coords.append([x, y])
        errors.append(error)

    coords = np.array(coords)
    errors = np.array(errors)

    # Triangulate element centers
    triangulation = tri.Triangulation(coords[:, 0], coords[:, 1])

    # Plot error using tripcolor
    plt.figure(figsize=(8, 6))
    tpc = plt.tripcolor(triangulation, errors, shading='gouraud', cmap=cmap)
    plt.colorbar(tpc, label="|Error|")
    plt.title(title)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis("equal")
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300)
        plt.close()
    else:
        plt.show()

# ------------------- RBF Mapping Main Function -------------------
def rbf_map_field(origin_model_part, destination_model_part, analytical_function, rbf_function, variable=Kratos.TEMPERATURE, support_radius=None):
    origin_coords = []
    origin_values = []
    destination_coords = []

    for element in origin_model_part.Elements:
        element_geometry = element.GetGeometry()
        x, y, z = element_geometry.Center().X, element_geometry.Center().Y, element_geometry.Center().Z
        value = analytical_function(x, y, z)
        origin_coords.append([x, y, z])
        origin_values.append(value)
        element.SetValue(variable, value)

    for element in destination_model_part.Elements:
        element_geometry = element.GetGeometry()
        x, y, z = element_geometry.Center().X, element_geometry.Center().Y, element_geometry.Center().Z
        destination_coords.append([x, y, z])

    origin_coords = np.array(origin_coords)
    origin_values = np.array(origin_values)
    R = support_radius or estimate_support_radius(origin_coords)

    A, n, poly_dim = build_origin_rbf_matrix(origin_coords, rbf_function, R)
    C = build_origin_destination_rbf_matrix(origin_coords, destination_coords, rbf_function, R)
    b = np.concatenate((origin_values, np.zeros(poly_dim + 1)))

    mapping_matrix = C @ np.linalg.inv(A)
    mapping_matrix_reduced = mapping_matrix[:, :n]
    print(mapping_matrix_reduced.shape)
    mapped_values = mapping_matrix_reduced @ b[:n]

    for i, element in enumerate(destination_model_part.Elements):
        element.SetValue(variable, mapped_values[i])

    compute_rmse(destination_model_part, analytical_function, variable)

    return mapping_matrix_reduced

# ------------------- RBF Mapping Function with transpose -------------------
def map_forces_transpose(model_part, variable, mapping_matrix, destination_forces):
    """
    Maps forces from destination back to origin using the transpose of the mapping matrix.

    Parameters:
    - mapping_matrix: NumPy array (shape: [n_dest, n_origin])
    - destination_forces: NumPy array (length: n_dest), i.e., forces on the destination

    Returns:
    - origin_forces: NumPy array (length: n_origin)
    """
    if mapping_matrix.shape[0] != len(destination_forces):
        raise ValueError("Mismatch: mapping_matrix rows must match number of destination forces.")

    origin_forces = mapping_matrix.T @ destination_forces

    for i, element in enumerate(model_part.Elements):
        element.SetValue(variable, [origin_forces[i], 0.0, 0.0])


# --------------------Main program ---------------------------------------------
if __name__ == "__main__":
    # ------------------------------ Create the origin and destination model parts -------------------------------
    # Create a model
    model = Kratos.Model()

    # Create an origin and destination model parts
    origin_model_part = model.CreateModelPart("OriginModelPart")
    destination_model_part = model.CreateModelPart("DestinationModelPart")

    # Route to the geometry files
    origin_cad_json_file = "origin_geometry.cad.json"
    destination_cad_json_file = "destination_geometry.cad.json"

    # Read the geometries
    Kratos.CadJsonInput(origin_cad_json_file).ReadModelPart(origin_model_part)
    origin_surface = origin_model_part.GetGeometry(2)
    Kratos.CadJsonInput(destination_cad_json_file).ReadModelPart(destination_model_part)
    destination_surface = destination_model_part.GetGeometry(2)

    # Create quadrature_point_geometries in the origin and destination
    origin_quadrature_point_geometries = Kratos.GeometriesVector()
    origin_surface.CreateQuadraturePointGeometries(origin_quadrature_point_geometries, 3)
    destination_quadrature_point_geometries = Kratos.GeometriesVector()
    destination_surface.CreateQuadraturePointGeometries(destination_quadrature_point_geometries, 3)

    element_id = 1
    shell_properties = origin_model_part.GetProperties()[1]
    shell_properties.SetValue(Kratos.THICKNESS, 0.1)
    shell_properties.SetValue(Kratos.YOUNG_MODULUS, 200000000)
    shell_properties.SetValue(Kratos.POISSON_RATIO, 0)
    shell_properties.SetValue(Kratos.CONSTITUTIVE_LAW, SMA.LinearElasticPlaneStress2DLaw())

    for i in range(0, len(origin_quadrature_point_geometries)):
        origin_model_part.CreateNewElement('Shell3pElement', element_id, origin_quadrature_point_geometries[i], shell_properties)
        element_id += 1

    element_id = 1
    for i in range(0, len(destination_quadrature_point_geometries)):
        destination_model_part.CreateNewElement('Shell3pElement', element_id, destination_quadrature_point_geometries[i], shell_properties)
        element_id += 1

    print(model)

    # Run RBF mapping
    mapping_matrix = rbf_map_field(
        origin_model_part,
        destination_model_part,
        analytical_function=analytical_field,
        rbf_function=wendland_c2_rbf,
        variable=Kratos.TEMPERATURE
    )
    
    # Define a force vector for the destination and map it to the origin
    destination_forces = []
    for element in destination_model_part.Elements:
        geom = element.GetGeometry()
        x, y, z = geom.Center().X, geom.Center().Y, geom.Center().Z
        destination_forces.append(analytical_force(x, y, z))
        element.SetValue(SMA.POINT_LOAD, [analytical_force(x, y, z), 0.0, 0.0])
    map_forces_transpose(origin_model_part, SMA.POINT_LOAD, mapping_matrix, np.array(destination_forces))

    # Plot the displacements and the error in the origin and destination
    tripcolor_plot(origin_model_part, Kratos.TEMPERATURE, title="Origin Field", filename="origin_field_disp")
    tripcolor_plot(destination_model_part, Kratos.TEMPERATURE, title="Destination Field", filename="destination_field_disp")
    tripcolor_error_plot(destination_model_part, analytical_field, Kratos.TEMPERATURE, filename="error_plot_disp")

    # Plot the mapped forces in the origin and destination
    tripcolor_plot(destination_model_part, SMA.POINT_LOAD, filename="destination_forces")
    tripcolor_plot(origin_model_part, SMA.POINT_LOAD, filename="origin_forces")

    # Check mapping matrix consistency
    row_sums = mapping_matrix.sum(axis=1)
    max_dev = np.max(np.abs(row_sums - 1.0))
    print(f"[INFO] Max deviation from row sum = 1: {max_dev:.2e}")

    # Check the conservation properies of the mapper
    origin_total = compute_total_force(origin_model_part, SMA.POINT_LOAD)
    destination_total = compute_total_force(destination_model_part, SMA.POINT_LOAD)

    print("[INFO] Total force in origin model part:", origin_total)
    print("[INFO] Total force in destination model part:", destination_total)

    difference = np.linalg.norm(origin_total - destination_total)
    print(f"[INFO] ||Origin - Destination|| = {difference:.4e}")
