clearvars; clc; close all;
config_figures;
this_dir = "initial_interactions/";

if ~exist(this_dir, 'dir')
    mkdir(this_dir)
end

breakcolor = "#A2142F";
coloru  = "#53A966";
colorux  = "#c3e166";
colorut = "#eaac8b";
colorerr = "#7E2F8E";
coloralt = "#0072BD";
coloraltx = "#4BDEEE";


% Specify the problem and its initial conditions. Epsilon must be symbolic
% to work with the perturbation expansion in nhyperbolic_characteristic_ICs
syms Sx Sepsilon Sux
A(Sux,Sepsilon) = Sux+Sepsilon*Sux^3/3;
epsilon = 0.5;
k=1.5;
u0(Sx) = sin(k*Sx)*exp(-Sx^2); u0t(Sx) = sym(0);

tic
n = 5; [L0,R0] = nhyperbolic_characteristic_ICs(A,n,u0,u0t);
disp("Derived initial conditions to order n=" + string(n+1))
% Note, we could further refine the initial conditions by using this as an
% initial guess with Newton's method
toc

Sc(Sux) = sqrt(diff(A(Sux,Sepsilon),Sux));
Q  = matlabFunction(subs(int(Sc,Sux),Sepsilon,epsilon));
c  = matlabFunction(subs(Sc,Sepsilon,epsilon));
cp = matlabFunction(subs(diff(Sc,Sux),Sepsilon,epsilon));
G = @(l,lx,r,rx) epsilon*(((r+2*l).*r.*lx + (2*r+l).*l.*rx) ./ (c(l)+c(r)));
L0 = subs(L0,Sepsilon,epsilon); R0 = subs(R0,Sepsilon,epsilon);

tfinal = 35; tswitch = 1.95; Nt = 10*floor(tfinal)+1; Nx = 80000; xbuf = 5;


% (equally spaced!) time steps split between MoL and transport_solver
t = linspace(0,tfinal,Nt)';
[~,i_split] = min(abs(t-tswitch));
tmol = t(1:i_split); ttrans = t(i_split:end);
dt = t(2)-t(1);

% Solve the potential system for r = R_x and l = L_x
par.PDE = @(h,t,u) potential_system(c,G,h,t,u);

par.ICs = {matlabFunction(diff(L0,argnames(L0))), matlabFunction(diff(R0,argnames(R0)))};
par.tspan  = tmol;
par.N      = Nx;
par.xbuf   = xbuf;
par.domain = [-par.tspan(end)-par.xbuf, par.tspan(end)+par.xbuf, par.N];
% the system of ODEs shouldn't be stiff
par.solver = @ode78; % ode78 appears to be the fastest for this problem
par.options = odeset(RelTol=1e-10, AbsTol=1e-10);

tic
[tmol,xmol,U] = MoL1D(par.solver, par.PDE, par.tspan, par.domain, par.ICs, par.options);
Lx1 = U(:,:,1); Rx1 = U(:,:,2);
disp("Finished computing with MoL")
toc
dx = xmol(2)-xmol(1);
L1 = cumtrapz(xmol,Lx1,2);
R1 = cumtrapz(xmol,Rx1,2);
ux1 = Lx1+Rx1;
u1  = L1+R1;
ut1 = gradient(u1,dt);

% transport_solver requires the initial conditions be functions (not lists
% of discrete points), so using the solution from MoL is a bit tricky. We
% will build Chebfuns from the equally spaced data and differentiate them

% TODO try pchip, makima, or spline to see if the random spikes disappear...
% interp_style = 'linear';
interp_style = 'pchip';
Rxinterp = griddedInterpolant(xmol,Rx1(end,:),interp_style,'nearest');
Rxxinterp = griddedInterpolant(xmol,gradient(Rx1(end,:),dx),interp_style,'nearest');


% Repeat for R
[Xr,tr,r2,S,tb] = transport_solver(Q,c,cp, ...
    @(x) Rxinterp(x),@(x)Rxxinterp(x),xmol,ttrans,false);
for b = tb'
    tbi = find(abs(tr-b) < 1e-14);
    if(~isempty(tbi)), Xr(tbi,:) = []; tr(tbi) = []; r2(tbi,:) = []; S(tbi,:) = []; end
end
R2 = NaN*ones(size(r2)); idxs = ~isnan(Xr(1,:)); Ijump = find(isnan(Xr(1,:)));
for i = 1:length(tr)
    R2(i,idxs) = cumtrapz(Xr(i,idxs),r2(i,idxs));
end
Rt2(:,idxs) = gradient(R2(:,idxs)', dt)' - gradient(Xr(:,idxs)', dt)'.*r2(:,idxs);
disp("Found long-time dynamics with transport_solver")



%%
f = figure; f.Position = [476 360 900 330];
xlimits = [-3.5,3.5];

tiledlayout(2,3)
nexttile([2,2])
nskip = floor(Nt/5);
waterfall(xmol, tmol(1), R1(1,:), LineWidth=2, FaceAlpha=0); hold on
waterfall(Xr(1:nskip:end,:)-tr(1:nskip:end),repmat(tr(1:nskip:end),1,size(Xr,2)),R2(1:nskip:end,:),LineWidth=2,HandleVisibility='off'); hold on
plot3(S(:,1)-tr,tr,R2(:,Ijump(1)-1),color=breakcolor)
plot3(S(:,2)-tr,tr,R2(:,Ijump(2)-1),color=breakcolor)
plot3(S(:,3)-tr,tr,R2(:,Ijump(3)-1),color=breakcolor)
hold off

xlim(xlimits); ylim tight; zlim tight;
xlabel("$\tilde x = x-t$"); ylabel("$t$");
legend("$R$","$S$",Location="northeast")
colormap summer
view(25,13)


nexttile
xuxmax = zeros(1,length(tr)); uxmax = xuxmax;
xuxmin = zeros(1,length(tr)); uxmin = xuxmin;
for i = 1:length(tr)
    [M,Iguy] = max(R2(i,:));
    xuxmax(i) = Xr(i,Iguy)-tr(i);
    uxmax(i) = M;
    [M,Iguy] = min(R2(i,:));
    xuxmin(i) = Xr(i,Iguy)-tr(i);
    uxmin(i) = M;
end
plot(xmol,R1(1,:),Color=[hex2rgb(coloru),0.6]); hold on
% plot(Xr(1,:)-tr(1),R2(1,:),Color=[hex2rgb(coloru),0.6]); hold on
pguy = plot(Xr(end,:)-tr(end) ,R2(end,:),Color=[hex2rgb(coloru),1]);
plot(Xr(end,Ijump-1)-tr(end), R2(end,Ijump-1),'x',color=breakcolor);
penv = plot(xuxmax,uxmax,xuxmin,uxmin,Color="k");

legend([pguy;penv(1)],["$R$","Max \& Min"],Location="northwest")
xticklabels({})
hold off
xlim(xlimits); ylim tight
title("Initial \& final profile of $R$ and $R_x$")

nexttile
xuxmax = zeros(1,length(t)); uxmax = xuxmax;
xuxmin = zeros(1,length(t)); uxmin = xuxmin;
for i = 1:i_split
    [M,Iguy] = max(Rx1(i,:));
    xuxmax(i) = xmol(Iguy)-tmol(i);
    uxmax(i) = M;
end
for i = 1:length(tr)
    [M,Iguy] = max(r2(i,:));
    xuxmax(i+i_split) = Xr(i,Iguy)-tr(i);
    uxmax(i+i_split) = M;
end
plot(xmol,Rx1(1,:),Color=[hex2rgb(colorux),0.6]); hold on
% plot(Xr(1,:)-tr(1),r2(1,:),Color=[hex2rgb(colorux),0.6]); hold on
pguy = plot(Xr(end,:)-tr(end) ,r2(end,:),Color=[hex2rgb(colorux),1]);
plot([Xr(end,Ijump-1);Xr(end,Ijump+1)]-tr(end), [r2(end,Ijump-1);r2(end,Ijump+1)],'--',color=[0.7,0.7,0.7]); 
penv = plot(xuxmax,uxmax,Color="k");

legend([pguy;penv],["$R_x$","Max"],Location="northwest")
xlabel("$\tilde x = x-t$")
hold off
xlim(xlimits); ylim tight
print(this_dir+'rwaterfall.eps','-depsc','-vector')



%% Plot profiles
f = figure; f.Position = [476 360 800 330];

% This is different bc hybrid_solver mangles the solution a bit with cell arrays...
tiledlayout(1,2)
nexttile
nskip = floor(Nt/2);
for i=[1,30,60]
    alpha = 0.4+0.6*i/Nt;
    plot(x{i},u{i},Color=[hex2rgb(coloru),alpha]); hold on
    % if i ~=1, plot(X(i,Ijump-1), U(i,Ijump-1),'x',color=breakcolor); end
end
hold off
% plot the breaks with interpolation
xlim([-4.5,4.5])

nexttile
for i=[1,30,60]
    alpha = 0.4+0.6*i/Nt;
    plot(x{i},ux{i},Color=[hex2rgb(colorux),alpha]); hold on
    % if i ~=1, plot(X(i,Ijump-1), U(i,Ijump-1),'x',color=breakcolor); end
end
hold off
xlim([-4.5,4.5])








% Ensure for u(x,t) that u(0,t) = u(L,t).
function [xfor,xint,xbac] = periodic_BCs(Nx)
    xfor = [2:Nx,1]; xint = 1:Nx; xbac = [Nx,1:Nx-1];
end

% Solve the system of PDE
%   l_t = F(l)l_x + G(l,l_x,r,r_x)
%   r_t = - F(r)r_x - G(l,l_x,r,r_x)
function u = potential_system(F,G,h,~,u)
% construct indices that ensure periodic boundary conditions are satisfied
u = reshape(u,[],2); Nx = size(u,1);
[xfor,xint,xbac] = periodic_BCs(Nx);

% Update the PDE on the interior of the domain
l = u(xint,1); lx = (u(xfor,1) - u(xbac,1))/(2*h);
r = u(xint,2); rx = (u(xfor,2) - u(xbac,2))/(2*h);

interaction = G(l,lx,r,rx);
u(xint,1) = F(l).*lx + interaction;
u(xint,2) = - F(r).*rx - interaction;

% reshape u_ut to build a long column vector as ode23 expects
u = reshape(u,Nx*2,1);
end



