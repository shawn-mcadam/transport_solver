function u = CL_conserving_implicit_stencil(epsilon,u0,u1,M,N,tau,h)
% Initialize the solution matrix
u = zeros(N, M);
u(:, 1) = u0;
u(:, 2) = u1;

% Neumann boundary conditions
u(1, :) = u(end-1, :); % because x(1) == x(N)
u(end, :) = u(2, :);

% Solver for the implicit scheme using fsolve
for n = 2:M-1
    u_prev = u(:, n-1);
    u_curr = u(:, n);

    % Define the nonlinear system of equations
    F = @(u_next) u_next - tau^2 * ((1 + epsilon*((u_next([2:end 2]) - u_next([end-1 1:end-1])) / (2*h)).^2) .* ...
                                    (u_next([2:end 2]) - 2*u_next + u_next([end-1 1:end-1])) / h^2) - 2*u_curr + u_prev;

    % Solve the nonlinear system using fsolve
    options = optimoptions('fsolve', 'Display', 'none', 'TolFun', 1e-6, 'TolX', 1e-6);
    u_next = fsolve(F, u_curr, options);

    % Enforce Neumann boundary conditions for u_next
    u(:, n+1) = u_next;
    u(1, n+1) = u(end-1, n+1);
    u(end, n+1) = u(2, n+1);
end

end