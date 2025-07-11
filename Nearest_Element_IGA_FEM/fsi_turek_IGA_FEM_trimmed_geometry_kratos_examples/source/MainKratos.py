import KratosMultiphysics as KM
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
        self.time_history = []
        self.disp_x_history_A = []
        self.disp_y_history_A = []

    def OutputSolutionStep(self):
        super().OutputSolutionStep()

        solver = self._GetSolver()

        if hasattr(solver, "model") and hasattr(solver.model, "solver_wrappers"):
            structure_solver = solver.model.solver_wrappers.get("structure")
            if structure_solver and hasattr(structure_solver, "model"):
                sub_model = structure_solver.model
                model_part_name = "IgaModelPart"

                if hasattr(sub_model, "HasModelPart") and sub_model.HasModelPart(model_part_name):
                    model_part = sub_model.GetModelPart(model_part_name)

        wet_interface_sub_model_part = model_part.GetSubModelPart("Load_3")
        current_time = wet_interface_sub_model_part.ProcessInfo[KratosMultiphysics.TIME]

        # Only do this the first time
        if not hasattr(self, "tracking_conditions_initialized"):
            self.tracking_conditions_initialized = True

            # Choose target locations
            self.point_A = (0.26, 0.2)  # Replace with your target coordinates

            def get_closest_condition(target_point):
                min_dist = float("inf")
                closest_cond = None
                for cond in wet_interface_sub_model_part.Conditions:
                    center = cond.GetGeometry().Center()
                    dist = (center.X - target_point[0])**2 + (center.Y - target_point[1])**2
                    if dist < min_dist:
                        min_dist = dist
                        closest_cond = cond
                return closest_cond

            self.condition_A = get_closest_condition(self.point_A)

        # --- Evaluate displacement at A ---
        condition = self.condition_A
        geom = condition.GetGeometry()
        N = geom.ShapeFunctionsValues()

        solution_x_A = 0.0
        solution_y_A = 0.0
        for i, node in enumerate(condition.GetNodes()):
            solution_x_A += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X) * N[0, i]
            solution_y_A += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y) * N[0, i]

        self.time_history.append(current_time)
        self.disp_x_history_A.append(solution_x_A)
        self.disp_y_history_A.append(solution_y_A)

    def Finalize(self):
        super().Finalize()

        import matplotlib.pyplot as plt

        # Enable LaTeX-style rendering
        plt.rcParams["text.usetex"] = True
        plt.rcParams["font.family"] = "serif"

        # Plot
        plt.figure()
        plt.plot(self.time_history, self.disp_x_history_A, label=r"$u_x$ at Point A")
        plt.plot(self.time_history, self.disp_y_history_A, label=r"$u_y$ at Point A")

        # Labels with LaTeX formatting
        plt.xlabel(r"\textbf{Time} [s]")
        plt.ylabel(r"\textbf{Displacement} $u_x$ [m]")
        plt.title(r"\textbf{Displacement vs Time (Nearest Neighbor Mapper)}", fontsize=14)

        # Grid and legend
        plt.grid(True)
        plt.legend(loc="best")

        # Save and show
        plt.savefig("displacement_vs_time.pdf", dpi=300, bbox_inches="tight")
        plt.show()

if __name__ == "__main__":

    with open("ProjectParametersCoSim.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    simulation = CustomCoSimulationAnalysis(parameters)
    simulation.Run()