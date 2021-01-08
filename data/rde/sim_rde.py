from __future__ import absolute_import
from clawpack import riemann
from clawpack.riemann.euler_with_efix_1D_constants import *
from clawpack import pyclaw

import numpy as np
import scipy.io as sio

def step_Euler(solver,state,dt):
    # ODE solved explicitly with 2-stage, 2nd-order Runge-Kutta method.
    dt2 = dt/2.
    q = state.q

    # Reaction rate parameters:
    h = state.problem_data['h'] # heat release
    uc = state.problem_data['uc'] # activation temperature
    tSink = state.problem_data['tSink'] # sink temperature for 'radiation'
    nu1 = state.problem_data['nu1'] # viscosity
    nu2 = state.problem_data['nu2'] # viscosity
    e = state.problem_data['e'] # emissivity or 'radiation' coefficient
    a = state.problem_data['a'] # activation energy
    s = state.problem_data['s'] # injector stiffness / refill time constant 
    ut = state.problem_data['ut']# "injection pressure"
    r = state.problem_data['r'] # injector activation steepness
    n = state.problem_data['n'] # loss index

    # Define gradient:
    Z = state.aux[0,:]
    dx = (state.grid.upper[0]-state.grid.lower[0])/state.grid.num_cells[0]
    
    # Simplified reaction rate (rr):
    rr = np.exp((q[0,:] - uc)/a)*(1.0 - Z)
    beta = s/(1.0+np.exp(r*(q[0,:] - ut)))
    loss = e*(tSink - q[0,:])*(q[0,:]**n)

    # Get 1st stage quantities:
    qstar = q[0,:] + dt2*((rr*h) + loss + nu1*uxx(q[0,:],dx))
    ZStar  = Z + dt2*(rr - Z*beta + nu2*uxx(Z,dx))

    # Update variables and gradient:
    rr = np.exp((qstar - uc)/a)*(1.0 - ZStar)
    beta = s/(1.0+np.exp(r*(qstar - ut)))
    loss = e*(tSink - qstar)*(qstar**n)

    q[0,:] = q[0,:] + dt*((rr*h) + loss + nu1*uxx(qstar,dx))
    state.aux[0,:] = Z + dt*(rr - ZStar*beta + nu2*uxx(ZStar,dx))
    state.aux[1,:] = ZStar*beta
    
def uxx(u,dx):
    # This function takes the second spatial derivative of the periodic argument using numpy's built-in gradient function.
    
    # Concatenate the array pre- and post- array:
    uuu = np.concatenate((u,u,u))
    
    # Compute the second spatial derivative:
    qx = np.gradient(uuu,dx)                     
    qxx = np.gradient(qx,dx)
    
    # We only want the middle portion:
    l = len(u)
    
    return qxx[l:2*l]



"""
Nucleation
"""
## Specify Riemann Solver and instantiate solver object:
solver = pyclaw.ClawSolver1D(riemann.burgers_1D_py.burgers_1D)
solver.limiters = pyclaw.limiters.tvd.vanleer
solver.kernel_language = 'Python'

# Set up boundary conditions and CFL:
solver.step_source = step_Euler       
solver.bc_lower[0] = pyclaw.BC.periodic
solver.bc_upper[0] = pyclaw.BC.periodic
solver.aux_bc_lower[0]=pyclaw.BC.extrap
solver.aux_bc_upper[0]=pyclaw.BC.extrap
solver.max_steps = 100000
solver.cfl_desired = 0.1

# Specify domain and fields, and initial conditions:
L = 2.0*np.pi
mx = 1000
x = pyclaw.Dimension(0.0,L,mx,name='x')
domain = pyclaw.Domain(x)
state = pyclaw.State(domain,num_eqn=1,num_aux=2)
xc = state.grid.x.centers
# Two sech pulses:
#state.q[0,:] = 1.0 + 1.5*((1.0/np.cosh(1.0*(xc-1.0)))**20.0 + (1.0/np.cosh(1.0*(xc-3.0)))**20.0)

# One sech pulse
state.q[0,:] = 1.5*((1.0/np.cosh(1.0*(xc-1.0)))**10.0)

state.problem_data['efix']=True
state.problem_data['L'] = L

# Half combustion IC:
state.aux[0,:] = 0.0*state.q[0,:] #+ 0.5

# zero energy flux to start tracker:
state.aux[1,:] = state.aux[0,:]*0.0

# Model constants:
state.problem_data['efix']=True
state.problem_data['L'] = L
state.problem_data['s'] = 3.5

#TW Sims:
state.problem_data['h'] = 1.0        # heat release, D_cj = 2*h
state.problem_data['ut'] = 0.5      # 1/2*'injection pressure' or cutoff for injection activation function
state.problem_data['uc'] = 1.1       # activation temperature
state.problem_data['tSink'] = 0.0    # sink value (keep at zero)
state.problem_data['nu1'] = 0.0   # chamber diffusivity
state.problem_data['nu2'] = 0.0   # plenum diffusivity
state.problem_data['e'] = 0.11       # emissivity or radiation coefficient
state.problem_data['a'] = 0.30       # activation energy
state.problem_data['r'] = 5.0        # 'Steepness' of activation function
state.problem_data['n'] = 1.0        # Loss linearity (0 = linear, 1 = quadratic)

# Set up PyClaw controller:
claw = pyclaw.Controller()
claw.tfinal = 100
claw.solution = pyclaw.Solution(state,domain)
claw.solver = solver
claw.num_output_times = 10000
claw.output_format = None
claw.keep_copy = True
claw.write_aux_always = True

# Run the simulation:
claw.run();


# store the data:
u = np.empty((mx,claw.num_output_times+1))
z = np.empty((mx,claw.num_output_times+1))
mdot = np.empty((mx,claw.num_output_times+1))

for j in range(len(claw.frames)):
    u[:,j] = claw.frames[j].q
    z[:,j] = claw.frames[j].aux[0,:]
    mdot[:,j] = claw.frames[j].aux[1,:]
    
t = np.linspace(0, claw.tfinal, claw.num_output_times+1)
dt = t[1]-t[0]

mask = np.nonzero(t>0)[0]
t = t[mask]

tau = t*2*state.problem_data['h']/L
u = u[:, mask]
z = z[:, mask]


data = {"u":u, "z":z, 'x':xc, 't':t,
        "D":2*state.problem_data['h'], "L":L,
        "s":state.problem_data['s'], 'nu':state.problem_data['nu1']}
sio.savemat('../nucleation.mat', data)



"""
Annihilation
"""
## Specify Riemann Solver and instantiate solver object:
solver = pyclaw.ClawSolver1D(riemann.burgers_1D_py.burgers_1D)
solver.limiters = pyclaw.limiters.tvd.vanleer
solver.kernel_language = 'Python'

# Set up boundary conditions and CFL:
solver.step_source = step_Euler       
solver.bc_lower[0] = pyclaw.BC.periodic
solver.bc_upper[0] = pyclaw.BC.periodic
solver.aux_bc_lower[0]=pyclaw.BC.extrap
solver.aux_bc_upper[0]=pyclaw.BC.extrap
solver.max_steps = 100000
solver.cfl_desired = 0.1


# Specify domain and fields, and initial conditions:
L = 2.0*np.pi
mx = 1000
x = pyclaw.Dimension(0.0,L,mx,name='x')
domain = pyclaw.Domain(x)
state = pyclaw.State(domain,num_eqn=1,num_aux=2)
xc = state.grid.x.centers
# Two sech pulses:
state.q[0,:] = 1.0 + 1.5*((1.0/np.cosh(1.0*(xc-1.0)))**10.0 + (1.0/np.cosh(1.0*(xc-3.0)))**10.0)

# One sech pulse
#state.q[0,:] = 1.5*((1.0/np.cosh(1.0*(xc-1.0)))**10.0)

state.problem_data['efix']=True
state.problem_data['L'] = L

# Half combustion IC:
state.aux[0,:] = 0.0*state.q[0,:] #+ 0.5

# zero energy flux to start tracker:
state.aux[1,:] = state.aux[0,:]*0.0

# Model constants:
state.problem_data['efix']=True
state.problem_data['L'] = L
state.problem_data['s'] = 1.0

#TW Sims:
state.problem_data['h'] = 1.0        # heat release, D_cj = 2*h
state.problem_data['ut'] = 0.5      # 1/2*'injection pressure' or cutoff for injection activation function
state.problem_data['uc'] = 1.1       # activation temperature
state.problem_data['tSink'] = 0.0    # sink value (keep at zero)
state.problem_data['nu1'] = 0.0   # chamber diffusivity
state.problem_data['nu2'] = 0.0   # plenum diffusivity
state.problem_data['e'] = 0.11       # emissivity or radiation coefficient
state.problem_data['a'] = 0.30       # activation energy
state.problem_data['r'] = 5.0        # 'Steepness' of activation function
state.problem_data['n'] = 1.0        # Loss linearity (0 = linear, 1 = quadratic)

# Set up PyClaw controller:
claw = pyclaw.Controller()
claw.tfinal = 150
claw.solution = pyclaw.Solution(state,domain)
claw.solver = solver
claw.num_output_times = 15000
claw.output_format = None
claw.keep_copy = True
claw.write_aux_always = True

# Run the simulation:
claw.run();


# store the data:
u = np.empty((mx,claw.num_output_times+1))
z = np.empty((mx,claw.num_output_times+1))
mdot = np.empty((mx,claw.num_output_times+1))

for j in range(len(claw.frames)):
    u[:,j] = claw.frames[j].q
    z[:,j] = claw.frames[j].aux[0,:]
    mdot[:,j] = claw.frames[j].aux[1,:]
    
t = np.linspace(0, claw.tfinal, claw.num_output_times+1)
dt = t[1]-t[0]

mask = np.nonzero(t>5)[0]
t = t[mask]

tau = t*2*state.problem_data['h']/L
u = u[:, mask]
z = z[:, mask]

data = {"u":u, "z":z, 'x':xc, 't':t,
        "D":2*state.problem_data['h'], "L":L,
        "s":state.problem_data['s'], 'nu':state.problem_data['nu1']}
sio.savemat('../annihilation.mat', data)
