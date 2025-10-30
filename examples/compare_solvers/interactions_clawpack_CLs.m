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

tfinal = 35; tswitch = 2; Nt = 25*floor(tfinal)+1; Nx = 80000; xbuf = 5;
[t,x,u,ux,ut,tb] = LR_sum_solver(Q,c,cp,G,L0,R0,tfinal,tswitch,Nt,Nx,xbuf);


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


%% How does this solution compare against Clawpack?
% Small time comparison of conservation laws compared to MoL solution
% The custom MoL solver requires u(0) and u(dt). We'll compute the u(dt)
% from the initial condition with a much finer mesh because otherwise we
% get a small wave moving left, skewing the amplitude of this solution

tic
data = pyrunfile("myelasticity.py", "output", problem=3,a=x{end}(1),b=x{end}(end),tf=tfinal,Ntsteps=int32(length(t)-1),ncells=750);
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
% figure
% for i=1:length(tclaw)
%     plot(xclaw, uclaw(i,:)); hold on
%     plot(x{i},u{i})
%     hold off
%     % xlim([0,Xr(end,end)])
%     xlim(tclaw(i)+[-2.5,2.5]); ylim tight
% 
%     legend("claw","hybrid")
%     title("$t="+tclaw(i)+"$")
%     pause(0.05)
% end







%%

xlim tight; ylim tight; zlim tight;
xlabel("$x$"); ylabel("$t$"); zlabel("Amplitude");
legend("$R$","$S$",Location="northeast")
colormap summer
view(25,13)


nexttile
xmax = zeros(1,length(t)); umax = xmax;
for i = 1:length(t)
    [M,Iguy] = max(U(i,:));
    xmax(i) = X(i,Iguy);
    umax(i) = M;
end
for i=floor(linspace(1,length(t),2))
    alpha = 0.4+0.6*i/length(t);
    pguy = plot(X(i,:),U(i,:),Color=[hex2rgb(coloru),alpha]); hold on
    if i ~=1, plot(X(i,Ijump-1), U(i,Ijump-1),'x',color=breakcolor); end
end
plotenvelope = plot(xmax,umax,Color="k");
legend([pguy,plotenvelope],["$R$","Envelope"])
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
    alpha = (i==1)*0.6+(i~=1);
    pguy = plot(X(i,:),Ux(i,:),Color=[hex2rgb(colorux),alpha]); hold on
    if i ~=1, plot([X(end,Ijump-1);X(end,Ijump+1)], [Ux(end,Ijump-1);Ux(end,Ijump+1)],'--',color=[0.7,0.7,0.7]); end
end
plotenvelope = plot(xuxmax,uxmax,xuxmin,uxmin,Color="k");
legend([pguy,plotenvelope(1)],["$R_x$","Envelope"])
xlabel("$x$")
hold off;
xlim tight; ylim tight

print(this_dir+'rwaterfall.eps','-depsc','-vector')



%% Conservation laws
fig = figure; fig.Position = [476 360 1100 340];
colclaw = "#0072BD"; colchar  = "#53A966";

% Compute the absolute error in the approximation compared to u_tt from Clawpack
rows = 4; cols = 3;
tiledguy = tiledlayout(rows,cols,"Padding","none","TileSpacing","tight");


% Amplitude:
nexttile(cols+1)
amp = max(abs(uclaw),[],2);
amp_burg = cellfun(@(ui) max(abs(ui),[],2),u);
plot(tclaw, amp, Color=colclaw); hold on; plot(t,amp_burg,Color=colchar); xline(tswitch,color='r'); xline(tb); hold off;
title("Amplitude: $\mathrm{max}_x u$")
xlim tight;
set(gca,'xtick',[])

% legend
nexttile(1)
plot(nan(5,1),Color=colclaw); hold on; plot(nan(5,1),Color=colchar); hold off; axis off;
lgd = legend("Clawpack","System split wrt characteristics","location","north","fontsize",14);


% absolute error between clawpack and R
nexttile(2*cols+1,[2,1])
abserr = zeros(size(t));
for i=1:length(abserr)
    id_spc = x{i} > 0;
    abserr(i) = max(abs(interp1(xclaw,uclaw(i,:),x{i}(id_spc))-u{i}(id_spc)));
end
plot(t,abserr,Color=colorerr); hold on; xline(tswitch,color='r'); xline(tb); hold off
title("Absolute error: $\mathrm{max}_x |u_{\mathrm{ref}}-u_{\mathrm{char}}|$")
legend("Error",Location="southeast")
xlim tight; ylim tight


% momentum
nexttile(2)
cl = arrayfun(@(i) trapz(xclaw, utclaw(i,:)), 1:length(tclaw));
plot(tclaw, cl - cl(1), Color=colclaw); hold on; xline(tb); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("Momentum: $\int u_t \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
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
title("Center of mass: $\int tu_t-u \,\mathrm{d}x" + sign + num2scistr(abs(cl(1))) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(3*cols+2)
cl_burg = arrayfun(@(i) trapz(x{i}, t(i).*ut{i}-u{i}),1:length(t));
plot(t, cl_burg - cl_burg(1),"Color",colchar); hold on; xline(tswitch,color='r'); xline(tb); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int tu_t-u \,\mathrm{d}x" + sign + num2scistr(abs(cl_burg(1))) + "$")
xlim tight; ylim tight


% total energy
nexttile(3)
cl = arrayfun(@(i) trapz(xclaw, (utclaw(i,:).^2 + uxclaw(i,:).^2)/2 + epsilon*uxclaw(i,:).^4/12),1:length(tclaw));
plot(tclaw, cl - cl(1),Color=colclaw); hold on; xline(tb); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("Energy: $\int (u_t^2+u_x^2)/2 + \epsilon u_x^4/12 \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(cols+3)
cl_burg = arrayfun(@(i) trapz(x{i}, ut{i}.^2/2 + ux{i}.^2/2 + epsilon*ux{i}.^4/12),1:length(t));
pburg=plot(t, cl_burg - cl_burg(1),Color=colchar); hold on; xline(tswitch,color='r'); xline(tb); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int (u_t^2+u_x^2)/2 + \epsilon u_x^4/12 \,\mathrm{d}x" + sign + abs(cl_burg(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])


% Unknown physical interpretation
nexttile(2*cols+3)
cl = trapz(xclaw, utclaw.*uxclaw, 2);
plot(tclaw, cl - cl(1),"Color",colclaw); hold on; xline(tb); hold off
if cl(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int u_tu_x \,\mathrm{d}x" + sign + abs(cl(1)) + "$")
xlim tight; ylim tight
set(gca,'xtick',[])

nexttile(3*cols+3)
cl_burg = arrayfun(@(i) trapz(x{i}, ut{i}.*ux{i}), 1:length(t));
plot(t, cl_burg - cl_burg(1),"Color",colchar); hold on; xline(tswitch,color='r'); xline(tb); hold off
if cl_burg(1) < 0, sign = " + "; else, sign = " - "; end
title("$\int u_tu_x \,\mathrm{d}x" + sign + abs(cl_burg(1)) + "$")
xlim tight; ylim tight

xlabel(tiledguy,"Time",interpreter='latex')
print(this_dir+'err_CLs.eps','-vector','-depsc');
