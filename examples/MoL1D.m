function [t,x,U_Ut] = MoL1D(solver, PDE, tspan, domain, ICs, options)
% solver = odeXX
% PDE = @(Nx,step,t,u_ut) ...
%   PDE must provide the discritization, boundary conditions, and PDE itself
% timespan given directly to the odesolver
% domain = [xmin, xmax, step]
% ICs = {@(x) ..., @(x) ..., ...} is a cell array of univariate functions.
%   its length determines the time order of the PDE

% Construct the domain
xmin = domain(1); xmax = domain(2); Nx = domain(3);
step = (xmax-xmin)/Nx; x = linspace(xmin,xmax,Nx);

% Apply initial values functions to the constructed domain
Torder = length(ICs);
conds = zeros([Nx,Torder]);
for i = 1:Torder
    fpi = ICs{i};
    conds(:,i) = fpi(x);
end

options = odeset(options, Vectorized="on");
[t,U_Ut] = solver(@(t,u_ut) PDE(step,t,u_ut), tspan, conds, options);

% Reshape u,ut,... into a 3D matrix
U_Ut = reshape(U_Ut, length(t), Nx, Torder);

end