# Kratos imports
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.CoSimulationApplication import CoSimIO

# External imports
import numpy as np
from colorama import Fore, Style
import matlab.engine

# Starting the Matlab engine and adding the path
eng = matlab.engine.start_matlab()
eng.addpath(eng.genpath(r'/home/camarotti/Software/bem_acoustic_solver_main'), nargout = 0)

def print_decorator(func):
    def wrapped(*args, **kwargs):
        print(f"Executing {func.__name__} ...")
        ret = func(*args, **kwargs)
        print(f"Finished executing {func.__name__}")
        return ret
    return wrapped

### CoSim Register
s_connection_name = ""
bem_acoustic_mesh_name = "acoustic_mesh"
exported_acoustic_nodal_forces = CoSimIO.DoubleVector()
imported_acoustic_displacements = CoSimIO.DoubleVector()
model = Kratos.Model()

# Initial time and time step
current_frequency = 10.0  
delta_frequency = 1.0 

def cosimio_check_equal(a, b):
    assert a == b

@print_decorator
def AdvanceInTime(info):
    """Advances time for the dynamic simulation using CoSimIO.Info()."""

    global current_frequency, delta_frequency  # Ensure time variables are global

    # Compute the new time step
    current_frequency += delta_frequency
    print(Fore.CYAN + f"frequency = {current_frequency}" + Style.RESET_ALL) 

    # Create CoSimIO.Info() settings
    settings = CoSimIO.Info()
    settings.SetString("identifier", "AdvanceInTime")
    settings.SetString("connection_name", s_connection_name)
    settings.SetDouble("current_time", current_frequency)  # Send updated time

    return CoSimIO.Info()  # Return a new CoSimIO.Info() object

@print_decorator
def InitializeSolutionStep(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "InitializeSolutionStep")
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@print_decorator
def Predict(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "Predict")
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@print_decorator
def SolveSolutionStep(info):
    # Inputs for the solution step in Matlab
    filename = 'rect_1x0.5.mat'
    knot_insert = matlab.double([1, 1])
    order_elevation = matlab.double([0, 0])

    # Convert to Python list
    acoustic_displacements_list = list(imported_acoustic_displacements)
    
    # Convert to MATLAB-compatible column vector (Nx1)
    acoustic_displacements_matlab = matlab.double(acoustic_displacements_list, size=(len(acoustic_displacements_list), 1))

    acoustic_force_vector, RHS_acoustic, LHS_Acoustic = eng.calcacousticCoSim(
        acoustic_displacements_matlab, current_frequency, filename, knot_insert, order_elevation, nargout = 3)
    ## Here we solve the acoustic domain
    # model_part = model.GetModelPart(dummy_fluid_mesh_name)

    # dummy_fluid_nodal_forces.clear()

    # num_nodes = len(model_part.Nodes)

    # # Call MATLAB function
    # flat_force_vector = eng.send_force(num_nodes)
    # print(flat_force_vector)

    # for node in model_part.Nodes:
    #     dummy_fluid_nodal_forces.append([0.0, 0.0, -0.01])
    
    # print("I am solving the fluid domain")
    return CoSimIO.Info()

@print_decorator
def FinalizeSolutionStep(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "FinalizeSolutionStep")
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@print_decorator
def OutputSolutionStep(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "OutputSolutionStep")
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@print_decorator
def ImportData(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", info.GetString("identifier"))
    
    # Here the displacements are imported from the structural solver
    return_info = CoSimIO.ImportData(settings, imported_acoustic_displacements)
    return CoSimIO.Info()

@print_decorator
def ExportData(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", info.GetString("identifier"))

    # Here the acoustic domain forces are exported to the structural solver
    acoustic_force_to_export = []

    for force in acoustic_nodal_forces:
        acoustic_force_to_export.extend(force)  # Append the actual force values

    data_to_be_send = CoSimIO.DoubleVector(fluid_force_to_export)
    return_info = CoSimIO.ExportData(settings, data_to_be_send)

    return CoSimIO.Info()

@print_decorator
def ImportMesh(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "info_for_test")
    settings.SetString("name_for_check", "ImportMesh")
    if (info.Has("identifier")):
        settings.SetString("identifier_control", info.GetString("identifier"))
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@print_decorator
def ExportMesh(info):
    # Exporting mesh to Kratos
    info = CoSimIO.Info()
    info.SetString("identifier", bem_acoustic_mesh_name)
    info.SetString("connection_name", s_connection_name)

    model_part = model.CreateModelPart(bem_acoustic_mesh_name)

    # Inputs for the preprocessing step in Matlab
    filename = 'rect_1x0.5.mat'
    knot_insert = matlab.double([1, 1])
    order_elevation = matlab.double([0, 0])

    degree, knot_vector, weights, control_points, element_conn = eng.preprocessingBEM(
        filename, knot_insert, order_elevation, nargout = 5)
    
    # Create nodes in the acoustic model part
    node_index = 1
    for i in range(len(control_points)):
        # Create nodes in the acoustic model part
        model_part.CreateNewNode(node_index, control_points[i][0], control_points[i][1], control_points[i][2])
        node_index += 1

  
    return_info = CoSimIO.ExportMesh(info, model_part)
    return info

# Connection Settings
settings = CoSimIO.Info()
settings.SetString("my_name", "call_acoustic_bem_solver")
settings.SetString("connect_to", "acoustic_bem_solver")
settings.SetInt("echo_level", 1)
settings.SetString("version", "1.25")
settings.SetString("communication_format", "file") 

# Connecting
return_info = CoSimIO.Connect(settings)
cosimio_check_equal(return_info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Connected)
s_connection_name = return_info.GetString("connection_name")

# registering the functions
fct_info = CoSimIO.Info()
fct_info.SetString("connection_name", s_connection_name)

fct_info.SetString("function_name", "AdvanceInTime")
CoSimIO.Register(fct_info,           AdvanceInTime)

fct_info.SetString("function_name", "InitializeSolutionStep")
CoSimIO.Register(fct_info,           InitializeSolutionStep)

fct_info.SetString("function_name", "Predict")
CoSimIO.Register(fct_info,           Predict)

fct_info.SetString("function_name", "SolveSolutionStep")
CoSimIO.Register(fct_info,           SolveSolutionStep)

fct_info.SetString("function_name", "FinalizeSolutionStep")
CoSimIO.Register(fct_info,           FinalizeSolutionStep)

fct_info.SetString("function_name", "OutputSolutionStep")
CoSimIO.Register(fct_info,           OutputSolutionStep)

fct_info.SetString("function_name", "ImportData")
CoSimIO.Register(fct_info,           ImportData)

fct_info.SetString("function_name", "ExportData")
CoSimIO.Register(fct_info,           ExportData)

fct_info.SetString("function_name", "ImportMesh")
CoSimIO.Register(fct_info,           ImportMesh)

fct_info.SetString("function_name", "ExportMesh")
CoSimIO.Register(fct_info,           ExportMesh)

# running the simulation
# externally orchestrated
run_info = CoSimIO.Info()
run_info.SetString("connection_name", s_connection_name)
CoSimIO.Run(run_info)

# Disconnecting
disconnect_settings = CoSimIO.Info()
disconnect_settings.SetString("connection_name", s_connection_name)
return_info = CoSimIO.Disconnect(disconnect_settings)
cosimio_check_equal(return_info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Disconnected)


# # Inputs for the preprocessing step in Matlab
# filename = 'rect_1x0.5.mat'
# knot_insert = matlab.double([1, 1])
# order_elevation = matlab.double([0, 0])

# degree, knot_vector, weights, control_points, element_conn = eng.preprocessingBEM(
#     filename, knot_insert, order_elevation, nargout = 5)
# print(degree)
# print(knot_vector)
# print(control_points)
# print(element_conn)
