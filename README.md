Method of characteristics solver for transport equations of the form $u_t - c'(u)u_x = 0$.

Some highlights:
- Computes breaks via the principle of equal areas. So, the breaks are introduced specifically to respect the conservation law $\frac{\mathrm{d}}{\mathrm{d}t}\int u \mathrm{d}t = 0$.
- Constructs the solution's mesh from the breaks: discontinuities are bracketed by two points separated by their machine epsilon.
- Breaks can merge

The example files test the accuracy in our computed solutions. This includes a comparison against a very accurate solver design specifically to respect conserved quanitites of a second-order equation $u_{tt}=(1+\epsilon u_x^2)u_{xx}$.
