%% Compare
% 1. Full MoL solution with conservation law preserving stencil
% 2. Burger's solver applied to the split system without the
% O(epsilon) interaction terms.
% The initial conditions are: $R(x,0)=e^{-x^2}, L(x,0) = 0$ so
% transport_solver is very close to exact
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
coloralt = "#0072BD";
coloraltx = "#4BDEEE";

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
tf = 35;
% tf = 3;


%% Push beyond the breaking time by using Clawpack
% The python script calls Clawpack. The corresponding Riemann solver is
% myelasticity_riemannsolver.py

tic
data = pyrunfile("myelasticity.py", "output", problem=1, tf=tf, Ntsteps=int64(3*tf), ncells=1000, a=-tf-6, b=tf+6);
toc
disp("Solved with Pyclaw.")
solution = data.frames.cell;

xclaw = solution{1}.state.grid.x.centers.double;
tclaw = cellfun(@(si) si.state.t,solution)';
solclaw = cell2mat(cellfun(@(si) si.state.q.double, solution,UniformOutput=false)');

uxclaw = solclaw(1:2:end,:);
utclaw = solclaw(2:2:end,:);
uclaw = cumtrapz(xclaw,uxclaw,2);

res = trapz(xclaw,uxclaw,2);



%% --------- High resolution, small time error/residual estimates ---------
% Using a stationary frame of reference for comparing MoL and MoC

Nsteps = 100000;
x0 = linspace(xclaw(1),xclaw(end),Nsteps);
t0 = tclaw;
dt = t0(2)-t0(1);

tic
[X,t,Ux,S,tbreaks,partition,XI] = transport_solver(Q,c,cp, ...
    matlabFunction(diff(R0,Sx)),matlabFunction(diff(R0,Sx,Sx)), ...
    x0,t0,false);

% Make it easier to compute ut and compare solutions by removing the two
% added tbreak timesteps
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
% Compute U_t as the directional derivative of U along characteristics
Ut(:,idxs) = gradient(U(:,idxs)', dt)' - gradient(X(:,idxs)', dt)'.*Ux(:,idxs);

% remove indices containing NaNs for numerics
xnumeric = X(:,idxs);
uxnumeric = Ux(:,idxs); unumeric = U(:,idxs); utnumeric = Ut(:,idxs);
toc
disp("Solved transport equation $r_t = \sqrt{r^2}r_x$")


% %%
% % Make a gif of the wave moving: u, u_x, and u_t
% % Also, test the clawpack soln actually looks like the characteristics soln
% f = figure; f.Position=[680 458 700 420];
% 
% pause(1e-5)
% mindelay = 0.015;
% [tsub,I] = subvector(t0,0.02); % NEVER forget this gif plays over tsub not t!
% gifname = 'uchar_vs_uclaw.gif';
% exportgraphics(f, gifname);
% for k = 1:length(tsub)
%     % plot(xnumeric(I(k),:),uxnumeric(I(k),:),color=colorux); hold on
%     % plot(xclaw,uxclaw(I(k),:),color=coloraltx); hold off
%     plot(xnumeric(I(k),:),unumeric(I(k),:),color=colorux); hold on
%     plot(xclaw,uclaw(I(k),:),color=coloraltx); hold off
% 
%     legend("Characteristic system", "Clawpack", Location="northwest")
%     % ylim tight; xlim(tsub(k)+[-4,4])
%     ylim tight; xlim tight
%     title("$t="+string(round(tsub(k),2))+"$")
% 
%     exportgraphics(f, gifname, Append=true);
% end
% movefile(gifname, this_dir);
% pause(1e-1)


%%
f = figure; f.Position=[680 458 900 280];
tiledlayout(2,2)
nexttile([2,1])
k=length(t0);
plot(X(end,:),Ux(k,:),color=coloru); hold on
for i=1:length(tbreaks)
    plot(X(end,[Ijump(i)-1,Ijump(i)+1])',Ux(k,[Ijump(i)-1,Ijump(i)+1])',"--",color=[0.7,0.7,0.7],HandleVisibility='off');
end
plot(xclaw,uxclaw(end,:),color=coloralt); hold off
ylim tight; xlim(t0(k)+[-2.5,3.2])
legend("Characteristic system", "Clawpack", Location="southwest")
title("$t="+t(end)+"$")

% The two zoom intervals have the same width: 0.37
nexttile
plot(X(end,:),Ux(k,:),color=coloru); hold on
for i=1:length(tbreaks)
    plot(X(end,[Ijump(i)-1,Ijump(i)+1])',Ux(k,[Ijump(i)-1,Ijump(i)+1])',"--",color=[0.7,0.7,0.7],HandleVisibility='off');
end
plot(xclaw,uxclaw(end,:),color=coloralt); hold off
xlim(t0(k)+[0.585,0.955])
title("Zoom into left break")

nexttile
plot(X(end,:),Ux(k,:),color=coloru); hold on
for i=1:length(tbreaks)
    plot(X(end,[Ijump(i)-1,Ijump(i)+1])',Ux(k,[Ijump(i)-1,Ijump(i)+1])',"--",color=[0.7,0.7,0.7],HandleVisibility='off');
end
plot(xclaw,uxclaw(end,:),color=coloralt); hold off
xlim(t0(k)+[2.63,3])
title("Zoom into right break")
print(this_dir+'clawpack_compare_breaks.eps','-vector','-depsc');


%% Compare conservation laws of the two solvers ----------
% Check if the two solvers respect some auxiliary conserved quantities of
% the PDE u_tt = A'(u_x)u_xx
fig = figure; fig.Position = [476 360 1000 340];
colclaw = "#0072BD"; colchar  = "#53A966";

% Compute the absolute error in the approximation compared to u_tt from Clawpack
rows = 4; cols = 3;
tiledguy = tiledlayout(rows,cols,"Padding","none","TileSpacing","tight");


% Amplitude:
nexttile(cols+1)
amp = max(abs(uclaw),[],2);
amp_burg = max(abs(unumeric),[],2);
plot(tclaw, amp, Color=colclaw); hold on; plot(t0,amp_burg,Color=colchar); xline(tbreaks); hold off;
title("Amplitude: $\mathrm{max}_x u$")
xlim tight;
set(gca,'xtick',[])

% legend
nexttile(1)
plot(nan(5,1),Color=colclaw); hold on; plot(nan(5,1),Color=colchar); hold off; axis off;
lgd = legend("Clawpack","System split wrt characteristics","location","north","fontsize",14);


% absolute error
% we must interpolate unumeric and uref to be over the same domain
nexttile(2*cols+1,[2,1])
abserr = zeros(size(t0));
% id0 = floor(length(xclaw)/2);
for i = 1:length(t0)
    abserr(i) = max(abs(interp1(xclaw,uclaw(i,:),xnumeric(i,:))-unumeric(i,:)));
    % abserr(i) = max(abs(interp1(xnumeric(i,:),unumeric(i,:),xclaw(id0:end))-uclaw(i,id0:end)));
end
plot(t0,abserr,Color=colorerr); hold on; xline(tbreaks); hold off
title("Absolute error: $\mathrm{max}_x|u_{\mathrm{ref}}-u_{\mathrm{char}}|$")
legend("Error",Location="southeast")
xlim tight; ylim tight


% momentum
nexttile(2)
cl = arrayfun(@(i) trapz(xclaw, utclaw(i,:)), 1:length(tclaw));
plot(tclaw, cl - cl(1), Color=colclaw); hold on; xline(tbreaks); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("Momentum: $\int u_t \,\mathrm{d}x" + sign + num2scistr(abs(cl(1))) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(cols+2)
cl_burg = arrayfun(@(i) trapz(xnumeric(i,:), utnumeric(i,:)), 1:length(t0));
plot(t0, cl_burg-cl_burg(1),Color=colchar); hold on; xline(tbreaks); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int u_t \,\mathrm{d}x" + sign + num2scistr(abs(cl_burg(1))) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])


% Center of mass
nexttile(2*cols+2)
cl = arrayfun(@(i) trapz(xclaw, tclaw(i)*utclaw(i,:)-uclaw(i,:)), 1:length(tclaw));
plot(tclaw, cl - cl(1),Color=colclaw); hold on; xline(tbreaks); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("Center of mass: $\int tu_t-u \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(3*cols+2)
cl_burg = arrayfun(@(i) trapz(xnumeric(i,:), t0(i).*utnumeric(i,:)-unumeric(i,:)),1:length(t0));
plot(t0, cl_burg - cl_burg(1),"Color",colchar); hold on; xline(tbreaks); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int tu_t-u \,\mathrm{d}x" + sign + abs(cl_burg(1)) + "$")
xlim tight; ylim tight


% total energy
nexttile(3)
cl = arrayfun(@(i) trapz(xclaw, (utclaw(i,:).^2 + uxclaw(i,:).^2)/2 + epsilon*uxclaw(i,:).^4/12),1:length(tclaw));
plot(tclaw, cl - cl(1),Color=colclaw); hold on; xline(tbreaks); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("Energy: $\int (u_t^2+u_x^2)/2 + \epsilon u_x^4/12 \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(cols+3)
cl_burg = arrayfun(@(i) trapz(xnumeric(i,:), utnumeric(i,:).^2/2 + uxnumeric(i,:).^2/2 + epsilon*uxnumeric(i,:).^4/12),1:length(t0));
pburg=plot(t0, cl_burg - cl_burg(1),Color=colchar); hold on; xline(tbreaks); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int (u_t^2+u_x^2)/2 + \epsilon u_x^4/12 \,\mathrm{d}x" + sign + abs(cl_burg(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])


% Unknown physical interpretation
nexttile(2*cols+3)
cl = trapz(xclaw, utclaw.*uxclaw, 2);
plot(tclaw, cl - cl(1),"Color",colclaw); hold on; xline(tbreaks); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int u_tu_x \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(3*cols+3)
cl_burg = arrayfun(@(i) trapz(xnumeric(i,:), utnumeric(i,:).*uxnumeric(i,:)),1:length(t0));
plot(t0, cl_burg - cl_burg(1),"Color",colchar); hold on; xline(tbreaks); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int u_tu_x \,\mathrm{d}x" + sign + abs(cl_burg(1)) + "$")
xlim tight; ylim tight

xlabel(tiledguy,"Time",interpreter='latex')
print(this_dir+'clawpack_err_CLs.eps','-vector','-depsc');


