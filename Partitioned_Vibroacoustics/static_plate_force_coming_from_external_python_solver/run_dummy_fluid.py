import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.CoSimulationApplication import CoSimIO
#import CoSimIO

import numpy as np


### CoSIM Register

s_connection_name = ""
dummy_fluid_mesh_name = "DUMMY_FLUID_MESH"
dummy_fluid_nodal_forces = []
#model_part = CoSimIO.ModelPart(dummy_fluid_mesh_name)
model = Kratos.Model()

def cosimio_check_equal(a, b):
    assert a == b

@time_decorator()
def AdvanceInTime(info):
    settings = CoSimIO.Info()
    settings.SetString("identifier", "AdvanceInTime")
    settings.SetString("connection_name", s_connection_name)
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def InitializeSolutionStep(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "InitializeSolutionStep")
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def Predict(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "Predict")
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def SolveSolutionStep(info):
    model_part = model.GetModelPart(dummy_fluid_mesh_name)

    dummy_fluid_nodal_forces.clear()

    for node in model_part.Nodes:
        dummy_fluid_nodal_forces.append([0.0, 0.0, -0.01])
    
    print("I am solving the fluid domain")
    return CoSimIO.Info()

@time_decorator()
def FinalizeSolutionStep(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "FinalizeSolutionStep")
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def OutputSolutionStep(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "OutputSolutionStep")
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def ImportData(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", info.GetString("identifier"))
    return CoSimIO.Info()

@time_decorator()
def ExportData(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", info.GetString("identifier"))

    fluid_force_to_export = []

    for force in dummy_fluid_nodal_forces:
        fluid_force_to_export.extend(force)  # Append the actual force values

    data_to_be_send = CoSimIO.DoubleVector(fluid_force_to_export)
    return_info = CoSimIO.ExportData(settings, data_to_be_send)

    return CoSimIO.Info()

@time_decorator()
def ImportMesh(info):
    settings = CoSimIO.Info()
    settings.SetString("connection_name", s_connection_name)
    settings.SetString("identifier", "info_for_test")
    settings.SetString("name_for_check", "ImportMesh")
    if (info.Has("identifier")):
        settings.SetString("identifier_control", info.GetString("identifier"))
    #CoSimIO.ExportInfo(settings)
    return CoSimIO.Info()

@time_decorator()
def ExportMesh(info):
    # Exporting mesh to Kratos
    info = CoSimIO.Info()
    info.SetString("identifier", dummy_fluid_mesh_name)
    info.SetString("connection_name", s_connection_name)

    model_part = model.CreateModelPart(dummy_fluid_mesh_name)

    # Create nodes 
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 0.5, 0.0, 0.0)
    model_part.CreateNewNode(3, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(4, 0.0, 0.5, 0.0)
    model_part.CreateNewNode(5, 0.5, 0.5, 0.0)
    model_part.CreateNewNode(6, 1.0, 0.5, 0.0)
    model_part.CreateNewNode(7, 0.0, 1.0, 0.0)
    model_part.CreateNewNode(8, 0.5, 1.0, 0.0)
    model_part.CreateNewNode(9, 1.0, 1.0, 0.0)

    # Create elements
    properties = model_part.CreateNewProperties(1)
    model_part.CreateNewElement("Element2D4N", 1, [1, 2, 5, 4], properties)
    model_part.CreateNewElement("Element2D4N", 2, [2, 3, 6, 5], properties)
    model_part.CreateNewElement("Element2D4N", 3, [4, 5, 8, 7], properties)
    model_part.CreateNewElement("Element2D4N", 4, [5, 6, 9, 8], properties)
  
    return_info = CoSimIO.ExportMesh(info, model_part)
    return info


# Connection Settings
settings = CoSimIO.Info()
settings.SetString("my_name", "run_dummy_fluid")
settings.SetString("connect_to", "dummy_fluid")
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