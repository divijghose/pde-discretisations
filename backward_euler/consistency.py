import numpy as np
import matplotlib.pyplot as pp
import os

def u_exact(x, t, a):
    return np.sin(2*np.pi*(x - a*t)) # An exact solution to PDE

def residual(x, t, dt, dx, a):
    # exact solution at n, n+1
    u_n   = np.sin(2*np.pi*(x - a*t))
    u_np1 = np.sin(2*np.pi*(x - a*(t + dt)))

    # exact solution at m+1, m-1
    u_np1_mp1 = np.roll(u_np1, -1)  
    u_np1_mm1 = np.roll(u_np1,  1)  

    # time derivative
    time_term = (u_np1 - u_n) / dt

    # centered spatial derivative
    space_term = a * (u_np1_mp1 - u_np1_mm1) / (2*dx)

    # Residual
    tau = np.abs(time_term + space_term)

    return tau

def consistency_dt():
    a = 1.0
    L = 1.0
    M = 400 
    dx = L / M # Fix dx small enough to have negligible effect on dt error
    x = np.arange(0, L, dx)
    t = 0.0

    dts = np.linspace(0.0001, 0.1, 100)
    errors = []

    for dt in dts:
        tau = residual(x, t, dt, dx, a)
        err = np.linalg.norm(tau, np.inf)
        errors.append(err)
    errors = np.array(errors)

    # Compute slope in log-log space 
    log_dt = np.log(dts)
    log_err = np.log(errors)
    slope, intercept = np.polyfit(log_dt, log_err, 1)

    results_dir = 'consistency_results'
    os.makedirs(results_dir, exist_ok=True)

    pp.figure(figsize=(6,4))
    pp.loglog(dts, errors, 'b-', linewidth=2, label='Residual') # Residuals 
    pp.loglog(dts, np.exp(intercept)*dts**slope, 'r--', label=f'Fit slope ≈ {slope:.2f}') # Reference line
    
    pp.xlabel(r'$\Delta t$')
    pp.ylabel('Residual (sup norm)')
    pp.legend()
    
    filepath = os.path.join(results_dir, "dt.png")
    pp.savefig(filepath)
    pp.close()

def consistency_dx():
    a = 1.0
    L = 1.0
    t = 0.0
    dt = 1e-7  # Fix dt small enough to have negligible effect on dx error

    Ms = np.arange(10, 1000) # Spatial refinements
    errors = []
    for M in Ms:
        dx = L / M
        x = np.linspace(0, L, M, endpoint=False)
        
        tau = residual(x, t, dt, dx, a)
        err = np.linalg.norm(tau, np.inf)
        errors.append(err)
    errors = np.array(errors)

    # Compute slope in log-log space
    log_dx = np.log(1.0 / Ms)  
    log_err = np.log(errors)
    slope, intercept = np.polyfit(log_dx, log_err, 1)

    results_dir = 'consistency_results'
    os.makedirs(results_dir, exist_ok=True)

    # Plot log-log
    pp.figure(figsize=(6,4))
    pp.loglog(L/Ms, errors, 'b-', linewidth=2, markersize=6, label='Residual')
    pp.loglog(L/Ms, np.exp(intercept)*(L/Ms)**slope, 'r--', label=f'Fit slope ≈ {slope:.2f}')
    
    pp.xlabel(r'$\Delta x$')
    pp.ylabel('Residual (sup norm)')
    pp.legend()
    
    filepath = os.path.join(results_dir, "dx.png")
    pp.savefig(filepath)
    pp.close()

consistency_dx()