import KratosMultiphysics as KM
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

"""
For user-scripting it is intended that a new class is derived
from CoSimulationAnalysis to do modifications
Check also "kratos/python_scripts/analysis_stage.py" for available methods that can be overridden
"""

parameter_file_name = "ProjectParametersA.json"
with open(parameter_file_name,'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())
    
model = KM.Model()
StructuralMechanicsAnalysis(model, parameters).Run()
