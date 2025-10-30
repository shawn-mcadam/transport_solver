%% Compare
% 1. Full MoL solution with conservation law preserving stencil
% 2. Burger's solver applied to the split system without the
% O(epsilon) interaction terms.
% The initial conditions are: $R(x,0)=e^{-x^2}, L(x,0) = 0$ so Burger's
% solver is very close to exact
clearvars; clc; close all;
config_figures;
this_dir = "interaction_free_dynamics/";
if ~exist(this_dir, 'dir')
    mkdir(this_dir)
end
breakcolor = "#A2142F";
coloru  = "#53A966";
colorux  = "#c3e166";
colorut = "#eaac8b";
colorerr = "#7E2F8E";
resolution = "-r1000";


% define the problem with symbolic variables
syms Sx Sepsilon Sux

% A(Sux,Sepsilon) = Sux + Sepsilon*Sux^3/3;
A(Sux,Sepsilon) = Sux + Sepsilon*Sux^3/3;
Sc(Sux) = sqrt(diff(A(Sux,Sepsilon),Sux));

% Derive the initial conditions
epsilon = 0.5;
R0(Sx)  = exp(-Sx^2);
L0(Sx)  = 0*Sx;

R0 = subs(R0,Sepsilon,epsilon);
r0 = matlabFunction(diff(R0,Sx));

% define Q,c,cp to be in a moving frame of reference
Q  = matlabFunction(subs(int(Sc,Sux),Sepsilon,epsilon));
c  = matlabFunction(subs(Sc,Sepsilon,epsilon));
cp = matlabFunction(subs(diff(Sc,Sux),Sepsilon,epsilon));

%% --------- High resolution, small time error/residual estimates ---------
% Using a stationary frame of reference for comparing MoL and MoC

Nsteps = 850;
tf = 4.05;
xbnd = 4;
x0 = linspace(-xbnd,xbnd,60*Nsteps);
t0 = linspace(0,tf,Nsteps)';
dt = t0(2)-t0(1);

tic
[X,t,Ux,S,tbreaks,partition,XI] = transport_solver(Q,c,cp, ...
    matlabFunction(diff(R0,Sx)),matlabFunction(diff(R0,Sx,Sx)), ...
    x0,t0,false);

% Make it easier to compute ut by removing the two added tbreak timesteps
for b = tbreaks'
    tbi = find(abs(t-b) < 1e-14);
    if(~isempty(tbi)), X(tbi,:) = []; t(tbi) = []; Ux(tbi,:) = []; S(tbi,:) = []; XI(tbi,:) = []; partition(tbi,:) = []; end
end
% [t,I] = subvector(t,dt);
% X = X(I,:); Ux = Ux(I,:); S = S(I,:); XI = XI(I,:); partition = partition(I,:);



Ijump = find(isnan(X(1,:))); U = NaN*ones(size(Ux)); Ut = NaN*ones(size(Ux));
idxs = ~isnan(X(1,:));
for i = 1:length(t)
    U(i,idxs) = cumtrapz(X(i,idxs),Ux(i,idxs));
end
% Compute U_t as the directional derivative of U
Ut(:,idxs) = gradient(U(:,idxs)', dt)' - gradient(X(:,idxs)', dt)'.*Ux(:,idxs);

% remove indices containing NaNs for numerics
xnumeric = X(:,idxs);
uxnumeric = Ux(:,idxs); unumeric = U(:,idxs); utnumeric = Ut(:,idxs);
toc
disp("Solved transport equation $r_t = \sqrt{r^2}r_x$ with a fine mesh")


% Small time comparison of conservation laws compared to MoL solution
% The custom MoL solver requires u(0) and u(dt). We'll compute the u(dt)
% from the initial condition with a much finer mesh because otherwise we
% get a small wave moving left, skewing the amplitude of this solution
N_IC = 32;
par.PDE = @(h,t,u_ut) quasilinear_2nd_order(@(ux) 1+epsilon*ux.^2,@periodic_BCs,h,t,u_ut);
par.ICs = {matlabFunction(R0), @(x) -Q(r0(x))};
% par.domain = [-xbnd-tf-1.2, xbnd+tf+1.2, N_IC*Nsteps];
par.domain = [-xbnd, xbnd+tf+1.2, N_IC*Nsteps];
par.solver = @ode23t; % (almost) unused in favour of CL_conserving_stensils
par.options = odeset(RelTol=3e-14, AbsTol=3e-14, Jacobian=jac_linear_sparsity(N_IC*Nsteps, 2));

Ntmol = 4; % Time refinement necessary for stability of MoL
tmol = linspace(0,t0(end),Ntmol*Nsteps)';
tau = tmol(2)-tmol(1);

tic
% Determine u1 by taking "one" time step with the method of lines (in
% quotes because Matlab's ODE solvers and hence MoL1D are adaptive in time)
[~,xmol,U] = MoL1D(par.solver, par.PDE, [0,tau], par.domain, par.ICs, par.options);
uref = U(:,:,1); u0 = uref(1,:); u1 = uref(end,:);

xmol = xmol(1:N_IC:end); u0 = u0(1:N_IC:end); u1 = u1(1:N_IC:end);
hmol = xmol(2) - xmol(1);
disp("Derived MoL initial conditions")
toc

tic
% Use a conservation law conserving stensil
M = length(tmol);
% uref    = CL_conserving_explicit_stencil(epsilon,u0,u1,M,length(xmol),tau,hmol);
uref = CL_conserving_semi_implicit_stencil(epsilon,u0,u1,M,length(xmol),tau,hmol);
% uref    = CL_conserving_implicit_stencil(epsilon,u0,u1,M,length(xmol),tau,hmol);
uref = permute(uref,[2 1]);
disp("Solved PDE with CL conserving MoL");
toc;

% compute finite difference approximations of u_x and u_t. Note that the
% conserved CLs are with respect to forward differences (not the O(h^2)
% centered differences). Also, ux^2 must be computed as uxref*uxref1.
% uxref = (uref(1:end-1,3:end)-uref(1:end-1,1:end-2))/(2*hmol);
uxref  = (uref(1:end-1,2:end)-uref(1:end-1,1:end-1))/hmol;
uxref1 = (uref(2:end,2:end)-uref(2:end,1:end-1))/hmol;
utref = (uref(2:end,:)-uref(1:end-1,:))/tau;

% Make a gif of the wave moving: u, u_x, and u_t
% Also, test that the MoL soln actually looks like the characteristics soln
%%
figure
pause(1e-5)
mindelay = 0.015;
[tsub,I] = subvector(t0,0.02); % NEVER forget this gif plays over tsub not t!
% gifname = 'smallt_R_Rx_Rt_profiles.gif';
gifname = 'uchar_vs_umol.gif';
for k = 1:length(tsub)
    % plot(X(I(k),:),Ux(I(k),:),color=colorux); hold on
    % plot(X(I(k),:),Ut(I(k),:),color=colorut);
    % plot(X(I(k),:), U(I(k),:),color=coloru);
    % for i = find(tbreaks < tsub(k))
    % plot(X(I(k),Ijump(i)-1), U(I(k),Ijump(i)-1),'x',color=breakcolor)
    % plot([X(I(k),Ijump(i)-1);X(I(k),Ijump(i)+1)], [Ux(I(k),Ijump(i)-1);Ux(I(k),Ijump(i)+1)],'--',color=[0.7,0.7,0.7])
    % plot([X(I(k),Ijump-1);X(I(k),Ijump+1)], [Ut(I(k),Ijump-1);Ut(I(k),Ijump+1)],'--',color=[0.7,0.7,0.7])
    % end
    % hold off
    % legend("$u_x$","$u_t$","$u$")

    plot(xnumeric(I(k),:),unumeric(I(k),:)); hold on
    plot(xmol,uref(Ntmol*I(k),:)); hold off
    legend("Approx. characteristic", "CL conserving MoL", Location="northwest")

    ylim tight; xlim tight
    title("$t="+string(round(tsub(k),2))+"$")
    if k == 1, gif(gifname,'nodither','DelayTime',mindelay); else, gif; end
end
movefile(gifname, this_dir);
pause(1e-1)

% Compare conservation laws of the two solvers ----------
% Check if the two solvers respect the conserved quantities of the PDE
% u_tt = A'(u_x)u_xx
% TODO Double check that the finite differences conservation laws are wrt
% centered space differences and 1st-order time differences

fig = figure; fig.Position = [476 360 1000 320];
colmol = "#0072BD"; colchar  = "#53A966";

% Compute the absolute error in the approximation compared to u_tt from MoL
% TODO pad the right side of this figures slightly...
r = 4; c = 3;
tiledguy = tiledlayout(r,c,"Padding","none","TileSpacing","tight");


% Amplitude:
nexttile(c+1)
amp = max(abs(uref),[],2);
amp_burg = max(abs(unumeric),[],2);
plot(tmol, amp, Color=colmol); hold on; plot(t0, amp_burg,Color=colchar); xline(tbreaks); hold off;
title("Amplitude: $\max_{x} u$")
xlim tight;
set(gca,'xtick',[])


% legend
nexttile(1)
plot(nan(5,1),Color=colmol); hold on; plot(nan(5,1),Color=colchar); hold off; axis off;
lgd = legend("CL converving MoL","System split wrt characteristics","location","north","fontsize",14);


% absolute error
% we must interpolate unumeric and uref to be over the same domain
nexttile(2*c+1,[2,1])
abserr = zeros(size(t0));
for i = 1:length(t0)
    abserr(i) = max(abs(interp1(xmol,uref(1+Ntmol*(i-1),:),xnumeric(i,:))-unumeric(i,:)));
end
plot(t0,abserr,Color=colorerr); hold on; xline(tbreaks); hold off
title("Absolute error: $\max_{x}|u_{\mathrm{ref}}-u_{\mathrm{char}}|$")
legend("Error",Location="southeast")
xlim tight; ylim tight


% momentum
nexttile(2)
cl = sum(hmol*utref, 2);
plot(tmol(1:end-1), cl - cl(1), Color=colmol); hold on; xline(tbreaks); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("Momentum: $\int u_t \,\mathrm{d}x" + sign + num2scistr(abs(cl(1))) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(c+2)
cl_burg = arrayfun(@(i) trapz(xnumeric(i,:), utnumeric(i,:)), 1:length(t0));
plot(t0, cl_burg-cl_burg(1),Color=colchar); hold on; xline(tbreaks); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int u_t \,\mathrm{d}x" + sign + num2scistr(abs(cl_burg(1))) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])


% Center of mass
nexttile(2*c+2)
cl = sum(hmol*(tmol(2:end).* utref - uref(2:end,:)), 2);
plot(tmol(1:end-1), cl - cl(1),Color=colmol); hold on; xline(tbreaks); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("Center of mass: $\int tu_t-u \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(3*c+2)
cl_burg = arrayfun(@(i) trapz(xnumeric(i,:), t0(i).*utnumeric(i,:)-unumeric(i,:)),1:length(t0));
plot(t0, cl_burg - cl_burg(1),"Color",colchar); hold on; xline(tbreaks); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int tu_t-u \,\mathrm{d}x" + sign + abs(cl_burg(1)) + "$")
xlim tight; ylim tight


% total energy
nexttile(3)
% cl = trapz(xmol(2:end-1), utref.^2/2 + uxref.^2/2 + epsilon*uxref.^4/12, 2);
cl = sum(hmol*(1/2*(uxref.*uxref1 + utref(:,2:end).^2)+ 1/12*epsilon*uxref.^2.*uxref1.^2), 2);
plot(tmol(1:end-1), cl - cl(1),Color=colmol); hold on; xline(tbreaks); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("Energy: $\int (u_t^2+u_x^2)/2 + \epsilon u_x^4/12 \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(c+3)
cl_burg = arrayfun(@(i) trapz(xnumeric(i,:), utnumeric(i,:).^2/2 + uxnumeric(i,:).^2/2 + epsilon*uxnumeric(i,:).^4/12),1:length(t0));
pburg=plot(t0, cl_burg - cl_burg(1),Color=colchar); hold on; xline(tbreaks); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int (u_t^2+u_x^2)/2 + \epsilon u_x^4/12 \,\mathrm{d}x" + sign + abs(cl_burg(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])


% Unknown physical interpretation
nexttile(2*c+3)
% cl = trapz(xmol(2:end-1), utref.*uxref, 2);
cl = sum(hmol*uxref.*utref(:,2:end), 2);
plot(tmol(1:end-1), cl - cl(1),"Color",colmol); hold on; xline(tbreaks); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("No known interpretation: $\int u_tu_x \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(3*c+3)
cl_burg = arrayfun(@(i) trapz(xnumeric(i,:), utnumeric(i,:).*uxnumeric(i,:)),1:length(t0));
plot(t0, cl_burg - cl_burg(1),"Color",colchar); hold on; xline(tbreaks); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int u_tu_x \,\mathrm{d}x" + sign + abs(cl_burg(1)) + "$")
xlim tight; ylim tight

xlabel(tiledguy,"Time",interpreter='latex')
print(this_dir+'err_CLs.eps','-vector','-depsc');

