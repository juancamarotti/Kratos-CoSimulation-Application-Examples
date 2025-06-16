import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend for remote execution

from geomdl import BSpline
from geomdl.visualization import VisMPL
import matplotlib.pyplot as plt

# Create a B-Spline curve instance
curve = BSpline.Curve()

# Set degree
curve.degree = 3

# Define control points
curve.ctrlpts = [[0, 0], [1, 2], [2, 2], [4, 0]]

# Define knot vector 
curve.knotvector = [0, 0, 0, 0, 1, 1, 1, 1]

# Evaluate the curve
curve.evaluate()

# Set up the visualization
curve.vis = VisMPL.VisCurve2D()

# Render the plot 
curve.render()

# Save the plot as a PNG file
plt.savefig("bspline_curve.png", dpi=300)
print("Plot saved as 'bspline_curve.png'")