Method of characteristics solver for transport equations of the form $u_t - c'(u)u_x = 0$.

Some highlights:
- Computes breaks via the principle of equal areas. So, each break is computed explicitly to respect the conservation law $\frac{\mathrm{d}}{\mathrm{d}t}\int u(x,t) \mathrm{d}t = 0.$
- Constructs the solution's mesh from the breaks: discontinuities are bracketed by two points separated by their machine epsilon.
- Breaks can merge.

The example files test the accuracy in our computed solutions. This includes a comparison against:
1. The implementation of the method of lines from "Invariant Conservation Law-Preserving Discretizations of Linear and Nonlinear Wave Equations" by Alexei et. al. This uses a stencil designed specifically to respect conserved quanitites of a second-order equation $u_{tt}=(1+\epsilon u_x^2)u_{xx}.$
2. TODO

![hippo](https://media4.giphy.com/media/v1.Y2lkPTc5MGI3NjExZXczdWNmM2V4OHp1YngwOGIyaXp1OHBnd3I0dndnZmZiOW9yMXZybiZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/DWDvrfAVFyJCcrsOC2/giphy.gif)
