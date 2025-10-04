#!/usr/bin/env python
# encoding: utf-8
r"""
Solitary wave formation in periodic nonlinear elastic media
===========================================================

Solve a one-dimensional nonlinear elasticity system:

.. math::
    \epsilon_t + u_x & = 0 \\
    (\rho(x) u)_t + \sigma(\epsilon,x)_x & = 0.

Here :math:`\epsilon` is the strain, :math:`\sigma` is the stress, 
u is the velocity, and :math:`\rho(x)` is the density.  
We take the stress-strain relation :math:`\sigma = e^{K(x)\epsilon}-1`;
:math:`K(x)` is the linearized bulk modulus.
Note that the density and bulk modulus may depend explicitly on x.


"""
from scipy.integrate import quad
import numpy as np

def qinit(state,c=1,M=0.5):
    x = state.grid.x.centers
    k = 1.5
    
    if problem==1: # exact travelling wave solution
        state.q[0,:] = -2.*x*np.exp(-x**2.)
        for i in range(len(state.q[0,:])):
            state.q[1,i] = -quad(lambda s: np.sqrt(c + M*s**2),0,state.q[0,i])[0]
    if problem==2: # O(epsilon) travelling wave solution
        state.q[0,:] = -2.*x*np.exp(-x**2.)
        state.q[1,:] = 2.*x*np.exp(-x**2.)
    if problem==3: # sine initial conditions
        #state.q[0,:] = k*(np.cos(k*x) - np.sign(x)*np.sin(k*x))*np.exp(-abs(x))
        state.q[0,:] = (k*np.cos(k*x)-2.*x*np.sin(k*x))*np.exp(-x**2.)
        state.q[1,:] = 0.


def setaux(x,rho=1,c=1,M=0.5):
    aux = np.empty([3,len(x)],order='F')
    # Density:
    aux[0,:] = rho
    # Bulk modulus:
    aux[1,:] = c
    # Nonlinearity
    aux[2,:] = M
    return aux


# solver_type=classic
def setup(use_petsc=0,solver_type='classic',outdir='./_output',
          tfinal=tf, num_output_times=Ntsteps):
    #from clawpack import riemann
    from myelasticity_riemannsolver import nonlinear_elasticity_1D

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    rs = nonlinear_elasticity_1D
    #rs = riemann.nonlinear_elasticity_1D_py.nonlinear_elasticity_1D

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)
        solver.char_decomp=0
    else:
        solver = pyclaw.ClawSolver1D(rs)

    solver.kernel_language = 'Python'

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic

    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap


    #xlower=a; xupper=b
    #xlower=-tfinal+2; xupper=tfinal+6.
    cells_per_layer=1000; mx=int(round(b-a))*cells_per_layer
    x = pyclaw.Dimension(a,b,mx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,solver.num_eqn,num_aux=3)

    #Set global parameters
    rho = 1.0
    KA  = 1.0
    KB  = 0.5

    #Initialize q and aux
    xc=state.grid.x.centers
    state.aux[:,:] = setaux(xc,rho,KA,KB)
    qinit(state)

    solver.max_steps = 5000000
    solver.fwave = True

    claw = pyclaw.Controller()
    claw.num_output_times = num_output_times
    claw.outdir = outdir
    claw.output_style = 1
    if outdir is None:
        claw.output_format = None
    claw.tfinal    = tfinal
    claw.solution  = pyclaw.Solution(state,domain)
    claw.solver    = solver
    claw.setplot   = setplot
    claw.keep_copy = True

    return claw

def setplot(plotdata):
    """
    I have no idea how to use Clawpack's plotting features and the docs are not
    particularly helpful...
    """
    plotdata.clearfigures() # Clear any old figures, axes, items data

    # Set plot_type to None to disable all plotting
    plotdata.plot_type = None 

    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    if 'problem' not in globals():
        problem = 1
    if 'tf' not in globals():
        tf = 35
    if 'a' not in globals():
        a = -tf-2
    if 'b' not in globals():
        b = tf+2
    if 'Ntsteps' not in globals():
        Ntsteps = int(tf)+1
    #if 'Nxsteps' not in globals():
    #    Nxsteps = 1000



    output = run_app_from_main(setup,setplot)
    #from clawpack.pyclaw import plot
    #plot.html_plot()
    #plot.interactive_plot()

