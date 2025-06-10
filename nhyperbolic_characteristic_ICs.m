function [L,R] = nhyperbolic_characteristic_ICs(A,n,u0,u0t)
% Determine L and R such that (L+R)(x,0) = u0(x) and (L+R)_t(x,0) = u0t(x)
% The result's spatial variable is the same as
% u0's and the perturbation variable is called Sepsilon. So
%   u(Sx,0,Sepsilon) \approx L(Sx,0,Sepsilon) + R(Sx,0,Sepsilon),
% where
%   L(Sx,0,Sepsilon) \approx L^0(Sx,0) + Sepsilon*L^1(Sx,0) + ...
%   R(Sx,0,Sepsilon) \approx R^0(Sx,0) + Sepsilon*R^1(Sx,0) + ...

% TODO could take c as the only argument instead of A

syms Sx Sepsilon R(Sx) L(Sx) s
syms r(Sx) [1,n]
assume(Sepsilon > 0)

% We need rvec so we can to refer to r1,r2,... with indices
% instead of doing function evaluation with the parentheses
rvec = r(Sx);

% set up the system of equations
cond1 = u0(Sx) == R(Sx) + L(Sx);
cond2 = u0t(Sx) == int(sqrt(diff(A(s,Sepsilon),s)),s,diff(R(Sx),Sx),diff(L(Sx),Sx));

% Sub cond1 into cond2, evaluate the integral, and Taylor expand R
cond2 = release(subs(cond2,L(Sx),rhs(isolate(cond1,L(Sx)))));
cond2 = subs(cond2,R(Sx),sum(r.*Sepsilon.^(0:n-1)));

% Solve for each Taylor coefficient. TODO differentiate RHS wrt Sepsilon!
for k=1:n
    orderk = limit(rhs(diff(cond2,Sepsilon,k-1)),Sepsilon,0,"right") == (k==1)*u0t(Sx);
    rk = dsolve(orderk, subs(rvec(k),Sx,-Inf) == 0);
    cond2 = subs(cond2,rvec(k),rk);
    rvec(k) = rk;
end
R(Sx) = sum(rvec.*Sepsilon.^(0:n-1));
L(Sx) = subs(rhs(isolate(cond1,L(Sx))));
end
