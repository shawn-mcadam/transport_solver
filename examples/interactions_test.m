clearvars; clc; close all;
config_figures;
this_dir = "initial_interactions/";

if ~exist(this_dir, 'dir')
    mkdir(this_dir)
end

% !! The following requires chebfun on MATLAB's path !!


syms Sx St Sepsilon Sux

A(Sux,Sepsilon) = Sux+Sepsilon*Sux^3/3;
u0(Sx) = exp(-Sx^2); u0t(Sx) = -diff(u0(Sx),Sx);
tic
n = 5; [L0,R0] = nhyperbolic_characteristic_ICs(A,n,u0,u0t);
disp("Derived initial conditions to order n=" + string(n+1))
toc

coloru = "#53A966";
colorux = "#c3e166";
colorut = "#eaac8b";
colorerr = "#7E2F8E";

tf = 2;
t0 = linspace(0,tf,100)';

% Set each parameter and give R & L to the numeric Burger's solver with
Sc(Sux) = sqrt(diff(A(Sux,Sepsilon),Sux));

epsilon = 0.5;
Q  = matlabFunction(subs(int(Sc,Sux),Sepsilon,epsilon));
c  = matlabFunction(subs(Sc,Sepsilon,epsilon));
cp = matlabFunction(subs(diff(Sc,Sux),Sepsilon,epsilon));

L0 = subs(L0,Sepsilon,epsilon); R0 = subs(R0,Sepsilon,epsilon);

% Order 1 equation directly. This is an integro-differential eqn when
% higher-order in epsilon
% F = @(ux) Q(ux);
% G = @(Lx,Rx) epsilon*Lx.*Rx.*(Rx+Lx)/2;

% Solve the potential system for r = R_x and l = L_x
G = @(l,lx,r,rx) epsilon*(((r+2*l).*r.*lx + (2*r+l).*l.*rx) ./ (c(l)+c(r)));

%%
% par.PDE = @(h,t,u) first_ord_system(Q,G,h,t,u);
% par.ICs = {matlabFunction(L0), matlabFunction(R0)};

par.PDE = @(h,t,u) potential_system(c,G,h,t,u);
par.ICs = {matlabFunction(diff(L0,Sx)), matlabFunction(diff(R0,Sx))};

par.tspan = t0;
par.N  = 30000;
par.xbuf = 4.25;
par.domain = [-par.tspan(end)-par.xbuf, par.tspan(end)+par.xbuf, par.N];
par.solver = @ode23;
par.options = odeset(RelTol=1e-10, AbsTol=1e-12);

tic
[tmol,xmol,U] = MoL1D(par.solver, par.PDE, par.tspan, par.domain, par.ICs, par.options);
Lx1 = U(:,:,1); Rx1 = U(:,:,2);
disp("Finished computing initial interactions with MoL")
toc
L1 = cumtrapz(xmol,Lx1,2);
R1 = cumtrapz(xmol,Rx1,2);


%%
[Xl,tl,l_approx,sl,tb_l] = transport_solver(@(ux) -Q(ux),@(ux) -c(ux),@(ux) -cp(ux), ...
    matlabFunction(diff(L0,Sx)),matlabFunction(diff(L0,Sx,Sx)), ...
    xmol,t0,false);
L_approx = NaN*ones(size(l_approx)); idxs = ~isnan(Xl(1,:));
for i = 1:length(tl)
    L_approx(i,idxs) = cumtrapz(Xl(i,idxs),l_approx(i,idxs));
end

[Xr,tr,r_approx,sr,tb_r] = transport_solver(Q,c,cp, ...
    matlabFunction(diff(R0,Sx)),matlabFunction(diff(R0,Sx,Sx)), ...
    xmol,t0,false);
R_approx = NaN*ones(size(r_approx)); idxs = ~isnan(Xr(1,:));
for i = 1:length(tr)
    R_approx(i,idxs) = cumtrapz(Xr(i,idxs),r_approx(i,idxs));
end


%%
f = figure; f.Position = [620 547 850 300];
tiledlayout(2,3)
nexttile([1,2])
plot(xmol,L1(1,:),Color=coloru); hold on
plot(xmol,L1(end,:),Color=colorux);
plot(Xl(end,:),L_approx(end,:),Color=colorut); hold off
ylabel("$L$"); xticklabels({})
legend("$L(x,0)$", "Full interations", "Interaction free after $t=0$")
xlim([-3.5,4]); ylim tight

nexttile([1,2])
plot(xmol,R1(1,:),Color=coloru); hold on
plot(xmol,R1(end,:),Color=colorux);
plot(Xr(end,:),R_approx(end,:),Color=colorut); hold off
ylabel("$R$"); xlabel("$x$")
legend("$R(x,0)$", "Full interations", "Interaction free after $t=0$",Location="northwest")
xlim([-3.5,4]); ylim tight


nexttile
errL = arrayfun(@(i) max(interp1(Xl(i,:),L_approx(i,:),xmol)-L1(i,:)),1:length(t0));
plot(t0,errL,Color=colorerr);
title("Max absolute error")
xticklabels({})
xlim tight; ylim tight;

nexttile
errR = arrayfun(@(i) max(interp1(Xr(i,:),R_approx(i,:),xmol)-R1(i,:)),1:length(t0));
plot(t0,errR,Color=colorerr);
xlabel("$t$")
xlim tight; ylim tight;

print(this_dir+'amp_error.eps','-vector','-depsc');



%% 

% Burgers_solver requires smooth initial conditions, so using the solution
% from MoL is a bit tricky. We will build Chebfuns from the equally spaced
% data and differentiate them
Lxcheb = chebfun(Lx1(end,:)',[xmol(1),xmol(end)],'equi');
Rxcheb = chebfun(Rx1(end,:)',[xmol(1),xmol(end)],'equi');

% how do the breaking times change for epsilon=0.5?
xi0 = linspace(par.domain(1),par.domain(2),300);
disp("Breaking time if $L$ was $0$")
[~,tb0] = burger_breakers(xi0,cp, matlabFunction(diff(u0)), matlabFunction(diff(u0,2)))
disp("Breaking time of $L$ with no interactions after $t=0$")
[~,tbL0] = burger_breakers(xi0,cp, matlabFunction(diff(L0)), matlabFunction(diff(L0,2)))
disp("Breaking time of $R$ with no interactions after $t=0$")
[~,tbR0] = burger_breakers(xi0,cp, matlabFunction(diff(R0)), matlabFunction(diff(R0,2)))
disp("Breaking time of $L$")
[~,tbL] = burger_breakers(xi0,cp, Lxcheb, diff(Lxcheb));
tbL = tbL+tf
disp("Breaking time of $R$")
[~,tbR] = burger_breakers(xi0,cp, Rxcheb, diff(Rxcheb));
tbR = tbR+tf



% Ensure for u(x,t) that u(0,t) = u(L,t).
function [xfor,xint,xbac] = periodic_BCs(Nx)
    xfor = [Nx,1:Nx-1]; xint = 1:Nx; xbac = [2:Nx,1];
end

% Solve the system of PDE:
%   L_t = - F(L_x) - G(L_x,R_x)
%   R_t = F(R_x) + G(L_x,R_x)
% Beyond O(epsilon^1), our system of PDEs of interest form an
% integro-differential equation and cannot be written like this
function u = first_ord_system(F,G,h,~,u)
% construct indices that ensure periodic boundary conditions are satisfied
u = reshape(u,[],2); Nx = size(u,1);
[xfor,xint,xbac] = periodic_BCs(Nx);
% Update the PDE on the interior of the domain
Lx = (u(xfor,1) - u(xbac,1))/(2*h);
Rx = (u(xfor,2) - u(xbac,2))/(2*h);
u(xint,1) = - F(Lx) - G(Lx, Rx);
u(xint,2) = F(Rx) + G(Lx, Rx);

% reshape u_ut to build a long column vector as ode23 expects
u = reshape(u,Nx*2,1);
end


% Solve the system of PDE
%   l_t = - F(l)l_x - G(l,l_x,r,r_x)
%   r_t =   F(r)r_x + G(l,l_x,r,r_x)
function u = potential_system(F,G,h,~,u)
% construct indices that ensure periodic boundary conditions are satisfied
u = reshape(u,[],2); Nx = size(u,1);
[xfor,xint,xbac] = periodic_BCs(Nx);

% Update the PDE on the interior of the domain
l = u(xint,1); lx = (u(xfor,1) - u(xbac,1))/(2*h);
r = u(xint,2); rx = (u(xfor,2) - u(xbac,2))/(2*h);

interaction = G(l,lx,r,rx);
u(xint,1) = - F(l).*lx - interaction;
u(xint,2) =   F(r).*rx + interaction;

% reshape u_ut to build a long column vector as ode23 expects
u = reshape(u,Nx*2,1);
end