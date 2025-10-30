function u = CL_conserving_semi_implicit_stencil(epsilon,u0,u1,M,N,tau,h)
% Initialize the solution matrix
u = zeros(N,M);
u(:, 1) = u0;
u(:, 2) = u1;

% Periodic boundary conditions
u(1, :) = u(end-1, :); % because x(1) == x(N)
u(end, :) = u(2, :);

% Solver for the implicit scheme using fsolve
for n = 2:M-1
    u_prev = u(:, n-1);
    u_curr = u(:, n);

    % calculate only for indices 2..end-1. Then 1<-end-1; end<-2
    % still do comp of N elemens, then impose periodicity
    uxx = (u_curr([2:end 1]) - 2*u_curr + u_curr([end 1:end-1])) / h^2;

    ux2 =  ( u_curr([2:end 1]) - u_curr ).^2/h^2;
    ux2m = ( u_curr - u_curr([end 1:end-1])).^2/h^2;

    Uprevx = ( u_prev([2:end 1]) - u_prev )/h;
    Uprevxm = ( u_prev - u_prev([end 1:end-1]) )/h;

    % objective function
    F = @(u_next) (u_next- 2*u_curr + u_prev)/tau^2 - uxx -1/6*...
        epsilon* ( ux2.*( (( u_next([2:end 1]) - u_next )/h) +Uprevx)-ux2m.*( (( u_next - u_next([end 1:end-1]) )/h)  +Uprevxm) )/h;

    % Solve the nonlinear system using fsolve
    options = optimoptions('fsolve', 'Display', 'none', 'TolFun', 1e-6, 'TolX', 1e-6);
    u_next = fsolve(F, u_curr, options);

    % Enforcing periodic boundary conditions for u_next
    u(:, n+1) = u_next;
    % u(1, n+1) = u(end-1, n+1);
    % u(end, n+1) = u(2, n+1);
end

end