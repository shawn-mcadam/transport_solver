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
% TODO could further refine this initial guess with Newton's method.
toc

Sc(Sux) = sqrt(diff(A(Sux,Sepsilon),Sux));
Q  = matlabFunction(subs(int(Sc,Sux),Sepsilon,epsilon));
c  = matlabFunction(subs(Sc,Sepsilon,epsilon));
cp = matlabFunction(subs(diff(Sc,Sux),Sepsilon,epsilon));
G = @(l,lx,r,rx) epsilon*(((r+2*l).*r.*lx + (2*r+l).*l.*rx) ./ (c(l)+c(r)));
L0 = subs(L0,Sepsilon,epsilon); R0 = subs(R0,Sepsilon,epsilon);

tfinal = 35; tswitch = 2; Nt = 25*floor(tfinal)+1; Nx = 30000; xbuf = 3.5;
[t,x,u,ux,ut,tb] = hybrid_solver(Q,c,cp,G,L0,R0,tfinal,tswitch,Nt,Nx,xbuf);



%% How does this solution compare against Clawpack?
% Small time comparison of conservation laws compared to MoL solution
% The custom MoL solver requires u(0) and u(dt). We'll compute the u(dt)
% from the initial condition with a much finer mesh because otherwise we
% get a small wave moving left, skewing the amplitude of this solution

tic
data = pyrunfile("myelasticity.py", "output", problem=3,a=x{end}(1),b=x{end}(end),tf=tfinal,Ntsteps=int32(length(t)-1));
% data = pyrunfile("myelasticity.py", "output", problem=3,a=-5,b=5,tf=3,Ntsteps=int8(25));
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


%%
figure
for i=1:length(tclaw)
    plot(xclaw, uclaw(i,:)); hold on
    plot(x{i},u{i})
    hold off
    % xlim([0,Xr(end,end)])
    xlim(tclaw(i)+[-2.5,2.5]); ylim tight
    
    legend("claw","hybrid")
    title("$t="+tclaw(i)+"$")
    pause(0.05)
end

%%
% compare abs error between clawpack and R
abserr = zeros(size(t));
for i=1:length(abserr)
    id_spc = x{i} > 0;
    abserr(i) = max(abs(interp1(xclaw,uclaw(i,:),x{i}(id_spc))-u{i}(id_spc)));
end
plot(t,abserr,color=colorerr)
title("Absolute error: $\max_{x}|u_{\mathrm{ref}}-u_{\mathrm{char}}|$")
legend("Error",Location="southeast")
xlim tight; ylim tight


%% Conservation laws








fig = figure; fig.Position = [476 360 1000 320];
colclaw = "#0072BD"; colchar  = "#53A966";

% Compute the absolute error in the approximation compared to u_tt from MoL
% TODO pad the right side of this figures slightly if possible...
rows = 4; cols = 3;
tiledguy = tiledlayout(rows,cols,"Padding","none","TileSpacing","tight");


% Amplitude:
nexttile(cols+1)
amp = max(abs(uclaw),[],2);
amp_burg = cellfun(@(ui) max(abs(ui),[],2),u);
plot(tclaw, amp, Color=colclaw); hold on; plot(t,amp_burg,Color=colchar); xline(tswitch,color='r'); xline(tb); hold off;
title("Amplitude: $\max_{x} u$")
xlim tight;
set(gca,'xtick',[])

% legend
nexttile(1)
plot(nan(5,1),Color=colclaw); hold on; plot(nan(5,1),Color=colchar); hold off; axis off;
lgd = legend("Clawpack","System split wrt characteristics","location","north","fontsize",14);


% absolute error
% we must interpolate unumeric and uref to be over the same domain
nexttile(2*cols+1,[2,1])


plot(t,abserr,Color=colorerr); hold on; xline(tswitch,color='r'); xline(tb); hold off
title("Absolute error: $\max_{x}|u_{\mathrm{ref}}-u_{\mathrm{char}}|$")
legend("Error",Location="southeast")
xlim tight; ylim tight


% momentum
nexttile(2)
cl = arrayfun(@(i) trapz(xclaw, utclaw(i,:)), 1:length(tclaw));
plot(tclaw, cl - cl(1), Color=colclaw);
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("Momentum: $\int u_t \,\mathrm{d}x" + sign + num2scistr(abs(cl(1))) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(cols+2)
cl_burg = arrayfun(@(i) trapz(x{i}, ut{i}), 1:length(t));
plot(t, cl_burg-cl_burg(1),Color=colchar); hold on; xline(tswitch,color='r'); xline(tb); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int u_t \,\mathrm{d}x" + sign + num2scistr(abs(cl_burg(1))) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])


% Center of mass
nexttile(2*cols+2)
cl = arrayfun(@(i) trapz(xclaw, tclaw(i)*utclaw(i,:)-uclaw(i,:)), 1:length(tclaw));
plot(tclaw, cl - cl(1),Color=colclaw); hold on; xline(tb); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("Center of mass: $\int tu_t-u \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(3*cols+2)
cl_burg = arrayfun(@(i) trapz(x{i}, t(i).*ut{i}-u{i}),1:length(t));
plot(t0, cl_burg - cl_burg(1),"Color",colchar); hold on; xline(tb); hold off
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
title("No known interpretation: $\int u_tu_x \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(3*cols+3)
cl_burg = arrayfun(@(i) trapz(xnumeric(i,:), utnumeric(i,:).*uxnumeric(i,:)),1:length(t0));
plot(t0, cl_burg - cl_burg(1),"Color",colchar); hold on; xline(tbreaks); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int u_tu_x \,\mathrm{d}x" + sign + abs(cl_burg(1)) + "$")
xlim tight; ylim tight

xlabel(tiledguy,"Time",interpreter='latex')
print(this_dir+'err_CLs.eps','-vector','-depsc');







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



