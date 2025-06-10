%% Compare
% 1. Full MoL solution with conservation law preserving stencil
% 2. Burger's solver applied to the split system without the
% O(epsilon) interaction terms.
% The initial conditions are: $R(x,0)=e^{-x^2}, L(x,0) = 0$ so Burger's
% solver is very close to exact
clearvars; clc; close all;
config_figures;
this_dir = "interaction_free_dynamics/";
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

% moving frame of reference
Q = @(ux) Q(ux) - ux;
c = @(ux) c(ux) - 1;


%% Set each parameter & solve for R
Nsteps = 80;
tf = 24.5;
x0 = linspace(-2,2.9,Nsteps);
t0 = linspace(0,tf,Nsteps)';

tic
[X,t,Ux,S,tbreaks,~,XI] = transport_solver(Q,c,cp, ...
    matlabFunction(diff(R0,Sx)),matlabFunction(diff(R0,Sx,Sx)), ...
    x0,t0,false);

Ijump = find(isnan(X(1,:))); U = NaN*ones(size(Ux));
idxs = ~isnan(X(1,:));
for i = 1:length(t)
    U(i,idxs) = cumtrapz(X(i,idxs),Ux(i,idxs));
end
toc
disp("Solved transport equation $r_t = \sqrt{r^2}r_x$")



%% Plot R with its breaks and some characteristics
% Get the char at the top of the domain & the ones bracketing each jump
f = figure; f.Position = [476 360 900 300];
tguy = tiledlayout(1,2,"Padding","none");

chars = XI(end,idxs) + c(r0(XI(end,idxs))).*t;
for i=2:length(t)
    chars(1:i,end+1:end+2*numel(Ijump)) = XI(i,[Ijump-1,Ijump+1]) + c(r0(XI(i,[Ijump-1,Ijump+1]))).*t(1:i);
    chars(i+1:end,end-2*numel(Ijump):end) = NaN;
end

nexttile
colorplot= pcolor(X,t,U,LineWidth=0.1); hold on
% surf(X,t,U,EdgeColor='none'); hold on
% plot(X,t,LineWidth=0.1,Color="k")
splot    = plot(S,t,"-",color=breakcolor,LineWidth=2); 
breakplt = yline(tbreaks); hold off
title("$R$'s computational mesh")
legend([colorplot,splot(1),breakplt(1)],["$R$","$S$","Breaks"],Location="northwest")

nexttile
pcolor(X,t,U,EdgeColor='none'); hold on
plot(chars,t,LineWidth=0.1,Color="k");
plot(S,t,"-",color=breakcolor,LineWidth=2);
yline(tbreaks);
title("$R$'s characteristics")
yticklabels({})

colormap summer
colorbar
xlim tight; ylim tight;
hold off

xlabel(tguy,"$x$",Interpreter="latex")
ylabel(tguy,"$t$",Interpreter="latex")

print(this_dir+'small_t_rchars.eps','-vector','-depsc')



%% large-time U and Ux
Nsteps = 300;
tf = 250;
x0 = linspace(-3,7,Nsteps/3);
t0 = linspace(0,tf,Nsteps)';

tic
[X,t,Ux,S,~,partition,XI] = transport_solver(Q,c,cp, ...
    matlabFunction(diff(R0,Sx)),matlabFunction(diff(R0,Sx,Sx)), ...
    x0,t0,false);

Ijump = find(isnan(X(1,:))); U = NaN*ones(size(Ux));
idxs = ~isnan(X(1,:));
for i = 1:length(t)
    U(i,idxs) = cumtrapz(X(i,idxs),Ux(i,idxs));
end
toc
disp("Solved large-time _low-res_ transport equation $r_t = \sqrt{r^2}r_x$")



%% Plot R with its breaks and some characteristics
% Get the char at the top of the domain & the ones bracketing each jump
chars = XI(end,idxs) + c(r0(XI(end,idxs))).*t;
for i=2:10:length(t)
    chars(1:i,end+1:end+2*numel(Ijump)) = XI(i,[Ijump-1,Ijump+1]) + c(r0(XI(i,[Ijump-1,Ijump+1]))).*t(1:i);
    chars(i+1:end,end-2*numel(Ijump):end) = NaN;
end

for j = 1:size(chars,2)
    k = sum(partition(1,:) <= chars(1,j));
    for i=1:length(t)
        if ~(partition(i,k) < chars(i,j) && chars(i,j) < partition(i,max(k+1,size(partition,2))))
            chars(i,j) = NaN;
        end
    end
end

f = figure; f.Position = [476 360 1000 300];
tguy = tiledlayout(1,2);
nexttile
pcolor(X,t,U,EdgeColor='none'); hold on
plot(chars,t,partition,t,LineWidth=0.2,Color="k");
plot(S,t,"-",color=breakcolor,LineWidth=2);
title("$R$'s characteristics")
ylabel("$t$")
xlim([-1.8,5.5]); ylim tight;
colormap summer; colorbar
hold off

nexttile
topchars = XI(end,idxs) + c(r0(XI(end,idxs))).*t;
plot(topchars,t,LineWidth=0.2,Color="k"); hold on
plot(S,t,LineWidth=2,Color=breakcolor); hold off
xlim([0,5.3]); ylim tight;
yticklabels({});
title("Characteristics must travel through a break")


xlabel(tguy,"$x$",Interpreter="latex")

print(this_dir+'large_t_rchars.eps','-vector','-depsc')
% print(this_dir+'large_t_rchars.jpg','-djpeg',resolution)


%% large-time U and Ux
Nsteps = 700;
tf = 250;
% Nsteps = 100;
% x0 = linspace(-3,7,3*Nsteps);
x0 = linspace(-3,7,70*Nsteps);
t0 = linspace(0,tf,Nsteps)';

tic
[X,t,Ux,S,tbreaks,partition,XI] = transport_solver(Q,c,cp, ...
    matlabFunction(diff(R0,Sx)),matlabFunction(diff(R0,Sx,Sx)), ...
    x0,t0,false);

Ijump = find(isnan(X(1,:))); U = NaN*ones(size(Ux)); idxs = ~isnan(X(1,:));
for i = 1:length(t)
    U(i,idxs) = cumtrapz(X(i,idxs),Ux(i,idxs));

    % % TODO use a chebfun to integrate over these chebyshev nodes. 
    % uxchebfun = chebfun({chebfun(Ux(1,1:Ijump(1)-2)',X(1,[1,Ijump(1)-2]), 'chebkind', 1),...
    % chebfun(Ux(1,Ijump(1)+1:Ijump(2)-2)',X(1,[Ijump(1)+1,Ijump(2)-2]), 'chebkind', 1),...
    % chebfun(Ux(1,Ijump(2)+1:end-2)',X(1,[Ijump(2)+1,end-2]), 'chebkind', 1)});
    % U(i,idxs) = arrayfun(@(x) integral(uxchebfun,X(i,1),x),X(i,idxs));
end
toc
disp("Solved large-time _high-res_ transport equation $r_t = \sqrt{r^2}r_x$")


%% waterfall plot shows the amplitude decays due to the break
f = figure; f.Position = [476 360 800 330];

tiledlayout(2,3)
nexttile([2,2])
nskip = floor(Nsteps/7);
% C = repmat(t(1:nskip:end),1,size(X,2));
waterfall(X(1:nskip:end,:),repmat(t(1:nskip:end),1,size(X,2)),U(1:nskip:end,:),LineWidth=2); hold on
plot3(S(:,1),t,U(:,Ijump(1)-1),color=breakcolor)
plot3(S(:,2),t,U(:,Ijump(2)-1),color=breakcolor)
hold off

xlim tight; ylim tight; zlim tight;
xlabel("$x$"); ylabel("$t$"); zlabel("Amplitude");
legend("$R$","$S$",Location="northeast")
% legend("Profiles","Envelope",Location="northeast")
colormap summer
% view(0,0)
view(18,13)


nexttile
xmax = zeros(1,length(t)); umax = xmax;
for i = 1:length(t)
    [M,Iguy] = max(U(i,:));
    xmax(i) = X(i,Iguy);
    umax(i) = M;
end
for i=floor(linspace(1,length(t0),2))
    alpha = (i==1)*0.6 + (i~=1);
    pguy = plot(X(i,:),U(i,:),Color=[hex2rgb(coloru),alpha]); hold on
    if i ~=1, plot(X(i,Ijump-1), U(i,Ijump-1),'x',color=breakcolor); end
end
plotenvelope = plot(xmax,umax,Color="k");
legend([pguy,plotenvelope],["R","Envelope"])
xticklabels({})
hold off

xlim tight; ylim tight
title("Initial \& final profile of $R$ and $R_x$")

nexttile
xuxmax = zeros(1,length(t)); uxmax = xuxmax;
xuxmin = zeros(1,length(t)); uxmin = xuxmin;
for i = 1:length(t)
    [M,Iguy] = max(Ux(i,:));
    xuxmax(i) = X(i,Iguy);
    uxmax(i) = M;
    [M,Iguy] = min(Ux(i,:));
    xuxmin(i) = X(i,Iguy);
    uxmin(i) = M;
end
for i=floor(linspace(1,length(t0),2))
    alpha = (i==1)*0.6 + (i~=1);
    pguy  = plot(X(i,:),Ux(i,:),Color=[hex2rgb(colorux),alpha]); hold on
    if i ~=1, plot([X(end,Ijump-1);X(end,Ijump+1)], [Ux(end,Ijump-1);Ux(end,Ijump+1)],'--',color=[0.7,0.7,0.7]); end
end
plotenvelope = plot(xuxmax,uxmax,xuxmin,uxmin,Color="k");
legend([pguy,plotenvelope(1)],["$R_x$","Envelope"])
xlabel("$x$")
hold off;
xlim tight; ylim tight

% print(this_dir+'rwaterfall.jpg','-djpeg',resolution)
print(this_dir+'rwaterfall.eps','-depsc','-vector')



%% Each break’s velocity compared to the characteristics that generate them
f = figure; f.Position = [476 360 850 280];
tguy = tiledlayout(2,5);
nexttile([1,3])
plot(t(2:end),diff(partition(:,2))./diff(t),color=breakcolor); hold on
plot(t,c(r0(XI(:,[Ijump(1)-1,Ijump(1)+1]))),"k"); hold off
legend("$S'(t)$","$c(r_0(\xi_+(t))$","$c(r_0(\xi_-(t))$")
title("Left break")
ylim tight; xlim tight
xticklabels({})
% fontsize(15,"points")

nexttile([1,3])
plot(t(2:end),diff(partition(:,3))./diff(t),color=breakcolor); hold on
plot(t,c(r0(XI(:,[Ijump(2)-1,Ijump(2)+1]))),"k"); hold off
title("Right break")
xlabel("$t$")
ylim tight; xlim tight

[~,general_entropy] = arrayfun(@(ul,ur) fminbnd(@(alpha) ...
    sign(ur-ul)*(Q(alpha*ur + (1-alpha)*ul) - alpha*Q(ur) - (1-alpha)*Q(ul)), ...
    0,1), r0(XI(:,Ijump(1)-1)), r0(XI(:,Ijump(1)+1)));
nexttile([2,2])
plot(t,general_entropy)
title("Generalized entropy condition: left break")
ylabel("$[q]$")
xlim([0,30]); ylim tight
xlabel("$t$")

ylabel(tguy,"Velocity",Interpreter="latex")

print(this_dir+'entropy.eps','-depsc','-vector')


% gen_ent = @(ul,ur,alpha) sign(ur-ul)*(Q(alpha*ur + (1-alpha)*ul) - alpha*Q(ur) - (1-alpha)*Q(ul));
% fplot(@(alpha) gen_ent(r0(XI(1,Ijump(1)-1)),r0(XI(1,Ijump(1)+1)),alpha),[0,1]);


%% Plot the residual error for the weak form of this PDE
fig = figure; fig.Position = [476 360 800 170];
tguy = tiledlayout(1,2);
nexttile
semilogx(t,U(:,end)-U(1,end)); hold on; xline(tbreaks); hold off
% title("Residual error for the conserved quantity $\displaystyle\int_{-\infty}^\infty r$")
title("$R(\infty,t) = \int_{-\infty}^\infty r(s,t)\, \mathrm{d}s$")
ylabel("Residual error")
xlabel("$t$")
xlim tight; ylim tight;

nexttile
semilogx(t,max(U,[],2)); hold on; xline(tbreaks); hold off
title("$\max_{x} R(x,t)$")
ylabel("Amplitude")
xlabel("$t$")
xlim tight; ylim tight;

print(this_dir+'residual_and_amplitude.eps','-depsc','-vector')


%% Make a gif of the wave moving:
figure
pause(1e-5)
mindelay = 0.025;
[tsub,I] = subvector(t0,0.25); % NEVER forget this gif plays over tsub not t!
gifname = 'R_Rx_profiles.gif';
for k = 1:length(tsub)
    plot(X(I(k),:), U(I(k),:),color=coloru); hold on
    plot(X(I(k),:),Ux(I(k),:),color=colorux);
    for i = find(tbreaks < tsub(k))
    plot(X(I(k),Ijump(i)-1), U(I(k),Ijump(i)-1),'x',color=breakcolor)
    plot([X(I(k),Ijump(i)-1);X(I(k),Ijump(i)+1)], [Ux(I(k),Ijump(i)-1);Ux(I(k),Ijump(i)+1)],'--',color=[0.7,0.7,0.7])
    end
    hold off

    ylim tight; xlim tight
    title("$t="+string(round(tsub(k),2))+"$")
    legend("$u$","$u_x$")

    if k == 1, gif(gifname,'nodither','DelayTime',mindelay); else, gif; end
end
movefile(gifname, this_dir);



%% --------- High resolution, small time error/residual estimates ---------
% Using a stationary frame of reference for comparing MoL and MoC
Q  = matlabFunction(subs(int(Sc,Sux),Sepsilon,epsilon));
c  = matlabFunction(subs(Sc,Sepsilon,epsilon));
cp = matlabFunction(subs(diff(Sc,Sux),Sepsilon,epsilon));

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



%% TODO move these functions to their own MoL related subdirectory & fix load_parameters
function S = jac_linear_sparsity(n,m)
    S = [sparse(zeros(n*(m-1),n)), speye(n*(m-1),n*(m-1));
         repmat(spdiags([1 -2 1],-1:1,n,n),1,m)];
end

% Ensure for u(x,t) that u(0,t) = u(L,t).
function [xfor,xint,xbac] = periodic_BCs(Nx)
    xfor = [Nx,1:Nx-1]; xint = 1:Nx; xbac = [2:Nx,1];
end

function u_ut = quasilinear_2nd_order(A_prime, boundaries, h, ~, u_ut1)
% odeXX gives u_ut1 as a long column vector. This code is easier to
% reason with if we reinterpret u_ut as a 2D matrix:
u_ut = reshape(u_ut1,[],2); Nx = size(u_ut,1);

% construct indices that ensure periodic boundary conditions are satisfied
[xfor,xint,xbac] = boundaries(Nx);
ux   = (u_ut(xfor,1) - u_ut(xbac,1))/(2*h);
uxx  = (u_ut(xfor,1) + u_ut(xbac,1) - 2*u_ut(xint,1))/h^2;

% Update the PDE on the interior of the domain
u_ut(xint,2) = A_prime(ux).*uxx;
% u_ut(xint,2) = (A(ux(xfor,1)) - A(ux(xbac,1)))/(2*h);

% copy ut into u's place
u_ut(:,1) = reshape(u_ut1(Nx+1:end),Nx,1);
% reshape u_ut to build a long column vector as ode23 expects
u_ut = reshape(u_ut,Nx*2,1);
end
