Method of characteristics solver for transport equations of the form $u_t - c'(u)u_x = 0$.

Some highlights:
- Computes breaks via the principle of equal areas. So, each break is computed explicitly to respect the conservation law $\frac{\mathrm{d}}{\mathrm{d}t}\int u(x,t) \mathrm{d}t = 0.$
- Constructs the solution's mesh from the breaks: discontinuities are bracketed by two points separated by their machine epsilon.
- Breaks can merge.

The example files test the accuracy in our computed solutions. This includes a comparison against:
1. The implementation of the method of lines from "Invariant Conservation Law-Preserving Discretizations of Linear and Nonlinear Wave Equations" by Alexei et. al. This uses a stencil designed specifically to respect conserved quanitites of a second-order equation $u_{tt}=(1+\epsilon u_x^2)u_{xx}.$
2. PyClaw's classical solver (second-order Godunov-type finite volume method based on Riemann solvers). The Riemann solver is adapted from [this example from Clawpack's gallary](https://www.clawpack.org/gallery/pyclaw/gallery/stegoton.html).


The following shows an example solution mesh compared to its characteristics. The solver carefully chooses mesh points so that discontinuities are represented with high accuracy.
![small_t_rchars-eps-converted-to](https://github.com/user-attachments/assets/b49fe44a-d822-4b41-9130-b0b0b842db1f)
