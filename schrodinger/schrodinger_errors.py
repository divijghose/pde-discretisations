from firedrake import *
import matplotlib.pyplot as plt

def dudt(u, u_, dt):
    return ((u - u_)/dt)

num_cells = 100 # Number of cells in each direction for the mesh
checkpointing = False # Checkpointing flag
t_max = 0.3
dt = 1/num_cells

k = 2 * pi
omega = k**2

# Create a unit square mesh
mesh = PeriodicIntervalMesh(num_cells, 1.0)
x = SpatialCoordinate(mesh)

V = VectorFunctionSpace(mesh, "CG", 1, dim=2)

v = TestFunction(V)

psi_0 = Function(V)
psi_0.interpolate(as_vector([cos(k*x[0]), sin(k*x[0])]))
psi_r0, psi_i0 = split(psi_0)

psi = Function(V)
psi.assign(psi_0)
psi_r, psi_i = split(psi)

F1 = (inner(dudt(psi_i, psi_i0, dt), v[0]) + inner(grad(psi_r), grad(v[0]))) * dx
F2 = (inner(dudt(psi_r, psi_r0, dt), v[1]) - inner(grad(psi_i), grad(v[1]))) * dx
F = F1 + F2

# step through time and record solution
t = 0.0
t_ufl = Constant(0.0)                 
psi_true = Function(V, name="psi_true")

times = []
l2_errors = []

while t < t_max - dt/2:
    t += dt
    solve(F == 0, psi)
    psi_0.assign(psi)


    t_ufl.assign(t)
    psi_true.interpolate(as_vector([cos(k * x[0] -  omega * t_ufl), sin(k * x[0] -  omega *  t_ufl)]))

    err_l2 = norms.errornorm(psi, psi_true)

    times.append(t)
    l2_errors.append(err_l2)

plt.plot(times, l2_errors)
plt.show()