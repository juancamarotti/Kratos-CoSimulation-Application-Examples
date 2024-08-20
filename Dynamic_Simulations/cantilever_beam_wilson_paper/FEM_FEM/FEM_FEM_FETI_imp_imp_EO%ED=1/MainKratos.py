import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

"""
For user-scripting it is intended that a new class is derived
from CoSimulationAnalysis to do modifications
Check also "kratos/python_scripts/analysis_stage.py" for available methods that can be overridden
"""

parameter_file_name = "cosim_fem_fem_small_2d_plate_feti_implicit_implicit_parameters.json"
with open(parameter_file_name,'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())
    

simulation = CoSimulationAnalysis(parameters)
simulation.Run()
