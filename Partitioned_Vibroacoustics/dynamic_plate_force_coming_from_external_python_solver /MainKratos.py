from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
import KratosMultiphysics.IgaApplication
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

# External imports
import importlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

import sys
import time

class CustomCoSimulationAnalysis(CoSimulationAnalysis):
    def __init__(self, parameters):
        """Call the base class constructor."""
        super().__init__(parameters)
    
    def OutputSolutionStep(self):
        super().OutputSolutionStep()

        solver = self._GetSolver()

        if hasattr(solver, "model") and hasattr(solver.model, "solver_wrappers"):
            structure_solver = solver.model.solver_wrappers.get("structure")  # Direct access

            if structure_solver and hasattr(structure_solver, "model"):
                sub_model = structure_solver.model
                model_part_name = "IgaBackgroundMeshModelPart"

                if hasattr(sub_model, "HasModelPart") and sub_model.HasModelPart(model_part_name):
                    model_part = sub_model.GetModelPart(model_part_name)
                
        # Lists to save the results
        x_coords = []
        y_coords = []
        solutions = []

        for element in model_part.Elements:
            element_geometry = element.GetGeometry()
                
            x = element_geometry.Center().X
            y = element_geometry.Center().Y

            index = 0
            solution = 0.0
            N = element_geometry.ShapeFunctionsValues()
            for node in element.GetNodes():
                nodal_solution = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)

                solution += nodal_solution * N[0,index]
                index += 1
                
            x_coords.append(x)
            y_coords.append(y)
            solutions.append(solution)

        # Convert lists to numpy arrays
        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)
        solutions = np.array(solutions)

        # Crear una cuadrícula para el contour plot
        grid_x, grid_y = np.mgrid[min(x_coords):max(x_coords):100j, min(y_coords):max(y_coords):100j]
        grid_solution = griddata((x_coords, y_coords), solutions, (grid_x, grid_y), method='linear')

        # Graficar el contour plot
        plt.figure(figsize=(8, 6))
        contour = plt.contourf(grid_x, grid_y, grid_solution, levels=20, cmap="viridis",vmin = min(solutions), vmax = max(solutions))
        plt.colorbar(contour, label="Solution (Temperature)")
        plt.xlabel("X Coordinate")
        plt.ylabel("Y Coordinate")
        plt.title("Contour Plot of Temperature Solution")

        # Guardar la figura en un archivo (por ejemplo, "contour_plot.png")
        plt.savefig("contour_plot.png", dpi=300, bbox_inches='tight', format='png')  # Cambia el formato si es necesario

        # Mostrar la figura (opcional)
        plt.show()

if __name__ == "__main__":

    with open("ProjectParametersCoSimMain.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    simulation = CustomCoSimulationAnalysis(parameters)
    simulation.Run()
