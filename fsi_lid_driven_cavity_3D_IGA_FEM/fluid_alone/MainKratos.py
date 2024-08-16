from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import sys
import time



if __name__ == "__main__":
    # Reading parameters from the *.json file
    with open("ProjectParametersCFD.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    # Creation of a model
    model = KratosMultiphysics.Model()

    # Generating the simulation
    simulation = FluidDynamicsAnalysis(model,parameters)

    # Running the simulation (exact sequence of steps can be seen by following the class hierarchy)
    simulation.Run()
