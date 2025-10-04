%% Compare
% 1. Full MoL solution with conservation law preserving stencil
% 2. Burger's solver applied to the split system without the
% O(epsilon) interaction terms.
% The initial conditions are: $R(x,0)=e^{-x^2}, L(x,0) = 0$ so Burger's
% solver is very close to exact
clearvars; clc; close all;
config_figures;
this_dir = "nonlinear_string_PDE/";
if exist(this_dir,'dir') ~= 7, mkdir(this_dir); end
breakcolor = "#A2142F";
coloru  = "#53A966";
colorux  = "#c3e166";
colorut = "#eaac8b";
colorerr = "#7E2F8E";
resolution = "-r1000";


% Define the problem with symbolic variables
syms Sx Sepsilon Sux

% Moving frame of reference removes the O(1) translation over time
A(Sux,Sepsilon) = atanh(Sux/sqrt(1+Sepsilon*Sux^2));
Sc(Sux) = simplify(sqrt(diff(A(Sux,Sepsilon),Sux)))-1; % the "-1" removes translation

% Derive the initial conditions
R0(Sx)  = exp(-Sx^2);
L0(Sx)  = 0*Sx;

% define Q,c,cp to be in a moving frame of reference
Q  = matlabFunction(int(Sc,Sux));
c  = matlabFunction(Sc);
cp = matlabFunction(diff(Sc,Sux));


%% Set each parameter & solve for R
Nsteps = 600;
tf = 535;
x0 = linspace(-6,6,Nsteps);
t0 = linspace(0,tf,Nsteps)';

f = figure; f.Position = [476 360 1000 310];
tguy = tiledlayout(1,3);

for epsilon=[0.5,0.73,0.85]
    tic
    [X,t,Ux,S,tbreaks,~,~] = transport_solver(@(ux) Q(ux,epsilon), ...
        @(ux) c(ux,epsilon), @(ux) cp(ux,epsilon), ...
        matlabFunction(diff(R0,Sx)),matlabFunction(diff(R0,Sx,Sx)), ...
        x0,t0,false);
    
    U = NaN*ones(size(Ux));
    for i = 1:length(t)
        idxs = ~isnan(X(i,:));
        U(i,idxs) = cumtrapz(X(i,idxs),Ux(i,idxs));
    end
    toc
    disp("Solved transport equation $r_t = \sqrt{r^2}r_x$")
    
    nexttile
    colorplot= pcolor(X,t,U,EdgeColor="none"); hold on
    splot    = plot(S,t,"-",color=breakcolor,LineWidth=2);
    breakplt = yline(tbreaks); hold off
    title("$\epsilon="+string(epsilon)+"$")
    if epsilon == 0.01
    legend([colorplot,splot(1),breakplt(1)],["$R$","$S$","Breaks"],Location="northwest")
    end
end

colormap summer; cb = colorbar; cb.Layout.Tile = 'east';
xlim tight; ylim tight;


xlabel(tguy,"$x$",Interpreter="latex")
ylabel(tguy,"$t$",Interpreter="latex")

print(this_dir+'spacetime_epsilons.eps','-depsc','-vector')



%% waterfall plot shows the amplitude decays due to the break
epsilon = 0.73;
x0 = linspace(-4,4,Nsteps);
t0 = linspace(0,300,Nsteps);
tic
[X,t,Ux,S,tbreaks,~,~] = transport_solver(@(ux) Q(ux,epsilon), ...
    @(ux) c(ux,epsilon), @(ux) cp(ux,epsilon), ...
    matlabFunction(diff(R0,Sx)),matlabFunction(diff(R0,Sx,Sx)), ...
    x0,t0,false);

Ijump = find(isnan(X(1,:))); U = NaN*ones(size(Ux));
for i = 1:length(t)
    idxs = ~isnan(X(i,:));
    U(i,idxs) = cumtrapz(X(i,idxs),Ux(i,idxs));
end
toc
disp("Solved transport equation $r_t = \sqrt{r^2}r_x$")


f = figure; f.Position = [476 360 800 330];

tiledlayout(2,3)
nexttile([2,2])
nskip = floor(Nsteps/7);
% C = repmat(t(1:nskip:end),1,size(X,2));
waterfall(X(1:nskip:end,:),repmat(t(1:nskip:end),1,size(X,2)),U(1:nskip:end,:),LineWidth=2); hold on
plot3(S(:,1),t,U(:,Ijump(1)-1),color=breakcolor)
plot3(S(:,2),t,U(:,Ijump(2)-1),color=breakcolor)
plot3(S(:,3),t,U(:,Ijump(3)-1),color=breakcolor)
plot3(S(:,4),t,U(:,Ijump(4)-1),color=breakcolor)
hold off

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


%% Further testing outside the range of currently supported parameters

% epsilon=0.689;
% 0.766
epsilon=0.72;
tic
r = matlabFunction(diff(R0,Sx));
rx = matlabFunction(diff(R0,Sx,Sx));
[X,t,Ux,S,tbreaks,partition,XI] = transport_solver(@(ux) Q(ux,epsilon), ...
    @(ux) c(ux,epsilon), @(ux) cp(ux,epsilon), ...
    r,rx, x0,t0,false);

U = NaN*ones(size(Ux));
for i = 1:length(t)
    idxs = ~isnan(X(i,:));
    U(i,idxs) = cumtrapz(X(i,idxs),Ux(i,idxs));
end
toc


figure
colorplot= pcolor(X,t,U,EdgeColor="none"); hold on
splot    = plot(partition,t,"-",color=breakcolor,LineWidth=2);
breakplt = yline(tbreaks); 
title("$\epsilon="+string(epsilon)+"$")
colormap summer

plot(XI(end,:)+c(r(XI(end,:)),epsilon).*t,t,"LineWidth",0.1,Color="k"); hold off


%%

gifname = 'epsilon_vs_spacetime.gif';
mindelay = 0.2;
pause(1e-5)
for epsilon = [0.65,0.66,0.675,0.685,0.69,0.7,0.71,0.72,0.73,0.743,0.75,0.755,0.768,0.77,0.774,0.7755,0.778,0.8]
    [X,t,Ux,S,tbreaks,~,~] = transport_solver(@(ux) Q(ux,epsilon), ...
        @(ux) c(ux,epsilon), @(ux) cp(ux,epsilon), ...
        matlabFunction(diff(R0,Sx)),matlabFunction(diff(R0,Sx,Sx)), ...
        x0,t0,false);
    U = NaN*ones(size(Ux));
    for i = 1:length(t)
        idxs = ~isnan(X(i,:));
        U(i,idxs) = cumtrapz(X(i,idxs),Ux(i,idxs));
    end
    
    pcolor(X,t,U,EdgeColor="none"); hold on
    plot(S,t,"-",color=breakcolor,LineWidth=2);
    yline(tbreaks); hold off
    title("$\epsilon="+string(epsilon)+"$")
    colormap summer
    ylim tight; xlim tight
    legend("$u(x,t)$","$S(t)$")

    if epsilon == 0.65, gif(gifname,'nodither','DelayTime',mindelay); else, gif; end
end
movefile(gifname, this_dir);


