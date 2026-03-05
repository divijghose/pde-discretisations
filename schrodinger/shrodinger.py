from firedrake import *
import matplotlib.pyplot as plt
# from firedrake.pyplot import tricontour


num_cells = 100 # Number of cells in each direction for the mesh
checkpointing = False # Checkpointing flag
# Create a unit square mesh
mesh = UnitSquareMesh(num_cells, num_cells)
if checkpointing:
    with CheckpointFile("poisson_example.h5", "w") as chk:
        chk.save_mesh(mesh)
        

# Define function space
V = FunctionSpace(mesh, "CG", 1)

# Define a test function and trial function
u = TrialFunction(V)
v = TestFunction(V)

# Define the manufactured solution and corresponding source term
u_analytical = Function(V)
x, y = SpatialCoordinate(mesh)
u_analytical.interpolate(-sin(omega * x) * sin(omega * y))
f = Function(V)
f.interpolate(-2 * (omega**2) * sin(omega * x) * sin(omega * y))
