function u = CL_conserving_explicit_stencil(epsilon,u0,u1,M,N,tau,h)
% Initialize the solution matrix
u = zeros(N, M);
u(:, 1) = u0;
u(:, 2) = u1;

% Periodic boundary conditions
u(N, :) = u(1, :);

% Solver for the implicit scheme using fzero
for n = 2:M-1
    u_prev = u(:, n-1);
    u_curr = u(:, n);
    
    for i = 2:N-1
        func = @(u_next_i) u_next_i - tau^2 * epsilon * (1 + ((u_curr(i+1) - u_curr(i-1)) / (2*h))^2) * (u_curr(i+1) - 2*u_curr(i) + u_curr(i-1)) / h^2 - 2*u_curr(i) + u_prev(i);
        u(i, n+1) = fzero(func, u_curr(i));
    end
    
    % Enforce periodic boundary conditions
    u(1, n+1) = u(N-1, n+1);
    u(N, n+1) = u(2, n+1);
end

end