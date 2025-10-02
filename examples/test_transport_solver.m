%% Test script for Burgers' solver
clearvars; clc; close all;
config_figures
is_test = "test_";
this_dir = is_test+"burgers/";
breakcolor = "#A2142F";

if exist(is_test+this_dir,'dir') ~= 7
    mkdir(is_test+this_dir)
end


%% --- Solve a common test for Burger's equation: u_t + u*u_x = 0 ---
% The flux
Q = @(u) u.^2/2; c = @(u) u; cp = @(u) 1;
% Initial condition and its derivative
phi = @(x) exp(-x.^2); phip = @(x) -2*x.*exp(-x.^2);


% Multivalued solution
Nmesh = 100; a = -6; b = 6;
x0 = linspace(a,b,Nmesh); t0 = linspace(0,10,Nmesh/5);
tic
[X,t,~,~,~] = transport_solver(Q,c,cp,phi,phip,x0,t0,true);
[~,tbreaks] = burger_breakers(x0,cp,phi,phip);
test1_time  = toc

figure
plot(X,t,color="black",LineWidth=0.5); hold on
yline(tbreaks);
xlim tight
hold off

% B&W sparsely populated characteristics of single valued solution
tic
[X,t,U,S,tbreaks,~,XI] = transport_solver(Q,c,cp,phi,phip,x0,t0,false);
test2_time = toc
disp("test 2 uses " + string(8*(numel(X)+numel(t)+numel(S)+numel(U))*1e-6) + "MB")

% compute the characteristics
idxs = ~isnan(X(1,:)); Ijump = find(isnan(X(1,:)));
chars = XI(end,idxs) + c(phi(XI(end,idxs))).*t;
for i=2:length(t)
    chars(1:i,end+1:end+2*numel(Ijump)) = XI(i,[Ijump-1,Ijump+1]) + c(phi(XI(i,[Ijump-1,Ijump+1]))).*t(1:i);
    chars(i+1:end,end-2*numel(Ijump):end) = NaN;
end

% plot the solution along with its characteristics
f = figure; f.Position = [680 458 800 380];
tiledlayout(1,2)
nexttile
plot(chars,t,color="black",LineWidth=0.5); hold on
plot(S,t,"x-",color=breakcolor,LineWidth=1,MarkerSize=1.75);
yline(tbreaks); hold off; xlim tight; ylim tight
nexttile
plot(X,t,color="black",LineWidth=0.5); hold on
plot(S,t,"x-",color=breakcolor,LineWidth=1,MarkerSize=1.75);
yline(tbreaks); hold off; xlim tight; ylim tight

% fine mesh of single valued solution
Nmesh = 350; a = -5; b = 7; tf = 15;
x0 = linspace(a,b,Nmesh); t0 = (linspace(0,tf,Nmesh)).^4/tf^3;
tic
[X,t,U,S,tbreaks] = transport_solver(Q,c,cp,phi,phip,x0,t0,false);
test3_time = toc
disp("test 3 uses " + string(8*(numel(X)+numel(U)+numel(t)+numel(S))*1e-6) + "MB")

figure
pcolor(X,t,U,EdgeColor='none'); hold on
plot(S,t,"-",color=breakcolor,LineWidth=2);
yline(tbreaks);
colormap summer
colorbar
xlim tight; ylim tight
hold off

figure
residual = arrayfun(@(i) trapz(X(i,~isnan(X(1,:))),U(i,~isnan(X(1,:)))),1:length(t));
plot(t,residual-residual(1)); hold on; xline(tbreaks); hold off
title("$\int_{-\infty}^\infty u(x,t) \mathrm{d}x - " + string(residual(1)) + "$");


%% tricky test example: u_t - 3u^2u_x = 0
% This weak solution quickly violates the entropy condition t~0.35
Q = @(u) -u.^3; c = @(u) -3*u.^2; cp = @(u) -6*u;
phi = @(x) -2*x.*exp(-x.^2); phip = @(x) 2*(2*x.^2-1).*exp(-x.^2);

Nmesh = 123; a = -4; b = 2;
x0 = linspace(a,b,Nmesh); t0 = linspace(0,5,Nmesh);
tic
[X,t,U,S,tbreaks] = transport_solver(Q,c,cp,phi,phip,x0,t0,false);
test4_time = toc
disp("test 4 uses " + string(8*(numel(X)+numel(U)+numel(t)+numel(S))*1e-6) + "MB")

figure
pcolor(X,t,U,EdgeColor='none'); hold on
plot(X,t,color="black",LineWidth=0.5);
plot(S,t,"-",color=breakcolor,LineWidth=2);
yline(tbreaks); hold off
colormap summer; colorbar; xlim tight; ylim tight

figure
surf(X,t,U,EdgeColor='none'); hold on
plot(S,t,"-",color=breakcolor,LineWidth=5); hold off
xlim tight; ylim tight
colormap summer;

figure
residual = arrayfun(@(i) trapz(X(i,~isnan(X(1,:))),U(i,~isnan(X(1,:)))),1:length(t));
plot(t,residual-residual(1)); hold on; xline(tbreaks); hold off
title("$\int_{-\infty}^\infty u(x,t) \mathrm{d}x - " + string(residual(1)) + "$");


%% I am curious: u_t - sin(u)u_x = 0
Q = @(u) -cos(u); c = @(u) sin(u); cp = @(u) cos(u);
% a=1.9; % only 1 break
% a=1.946; % barely two breaks
a=2.2; % difficult confluence
% a = 5; % Impossible break
phi = @(x) a*exp(-x.^2); phip = @(x) -2*a*x.*exp(-x.^2);

Nmesh = 102; a = -2.5; b = 6;
x0 = linspace(a,b,Nmesh);
t0 = linspace(0,9,Nmesh);
tic
[X,t,U,S,tbreaks,partition,XI] = transport_solver(Q,c,cp,phi,phip,x0,t0,false);
test5_time = toc
disp("test 4 uses " + string(8*(numel(X)+numel(U)+numel(t)+numel(S))*1e-6) + "MB")


% Plot R with its breaks and some characteristics
% Get the char at the top of the domain & the ones bracketing each jump
idxs = ~isnan(X(1,:)); Ijump = find(isnan(X(1,:)));
chars = XI(end,idxs) + c(phi(XI(end,idxs))).*t;
for i=2:length(t)
    chars(1:i,end+1:end+2*numel(Ijump)) = XI(i,[Ijump-1,Ijump+1]) + c(phi(XI(i,[Ijump-1,Ijump+1]))).*t(1:i);
    chars(i+1:end,end-2*numel(Ijump):end) = NaN;
end

figure
pcolor(X,t,U,EdgeColor='none'); hold on
plot(chars,t,LineWidth=0.1,Color="k");
plot(partition,t,"--",LineWidth=1,color=[0.6,0.6,0.6]);
plot(S,t,"-",color=breakcolor,LineWidth=2);
yline(tbreaks);
title("Characteristics")

colormap summer
colorbar
xlim tight; ylim tight;
hold off

figure
residual = arrayfun(@(i) trapz(X(i,~isnan(X(1,:))),U(i,~isnan(X(1,:)))),1:length(t));
plot(t,residual-residual(1)); hold on; xline(tbreaks); hold off
title("$\int_{-\infty}^\infty u(x,t) \mathrm{d}x - " + string(residual(1)) + "$");

%%
figure
pause(1e-5)
mindelay = 0.025;
gifname = 'exp_profile_sin_flux.gif';
for k = 1:length(t)
    plot(X(k,:), U(k,:)); hold on
    for i = find(tbreaks < t(k))
    plot([X(k,Ijump(i)-1);X(k,Ijump(i)+1)], [U(k,Ijump(i)-1);U(k,Ijump(i)+1)],'--',color=[0.7,0.7,0.7])
    end
    hold off

    xlim tight; ylim tight
    title("$t="+string(round(t(k),2))+"$")
    legend("$u$","Location","northwest")

    if k == 1, gif(gifname,'nodither','DelayTime',mindelay); else, gif; end
end
movefile(gifname, this_dir);



%% Periodic function with several breaks at the same times
Q = @(u) u.^2/2; c = @(u) u; cp = @(u) 1;
phi = @(x) sin(x); phip = @(x) cos(x);

Nmesh = 100; a = -2*pi; b = 2*pi;
x0 = linspace(a,b,Nmesh);
t0 = linspace(0,5,Nmesh);
tic
[X,t,U,S,tbreaks,~,XI] = transport_solver(Q,c,cp,phi,phip,x0,t0,false);
test4 = toc
disp("test 4 uses " + string(8*(numel(X)+numel(U)+numel(t)+numel(S))*1e-6) + "MB")


idxs = ~isnan(X(1,:)); Ijump = find(isnan(X(1,:)));
chars = XI(end,idxs) + c(phi(XI(end,idxs))).*t;
for i=2:length(t)
    chars(1:i,end+1:end+2*numel(Ijump)) = XI(i,[Ijump-1,Ijump+1]) + c(phi(XI(i,[Ijump-1,Ijump+1]))).*t(1:i);
    chars(i+1:end,end-2*numel(Ijump):end) = NaN;
end

f = figure;
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



%%
figure
pause(1e-5)
mindelay = 0.025;
gifname = 'sin_profile.gif';
for k = 1:length(t)
    plot(X(k,:), U(k,:)); hold on
    for i = find(tbreaks < t(k))
    plot([X(k,Ijump(i)-1);X(k,Ijump(i)+1)], [U(k,Ijump(i)-1);U(k,Ijump(i)+1)],'--',color=[0.7,0.7,0.7])
    end
    hold off

    xlim tight; ylim([-1,1])
    title("$t="+string(round(t(k),2))+"$")
    legend("$u$")

    if k == 1, gif(gifname,'nodither','DelayTime',mindelay); else, gif; end
end
movefile(gifname, this_dir);




%% Test solving uncoupled systems over compatible meshes.
%   u_{tt} = (A(u_x))_x
% is approximately the sum L + R where L and R move along the PDE's
% characteristics and do not interact. The initial conditions are
% constructed such that $u(x,0)=e^{-x^2}$ and $u_t(x,0)=-2u_x(x,0)$ to
% order epsilon
clearvars; clc

epsilon = 0.5;
Q  = @(ux) (asinh(sqrt(epsilon)*ux)/sqrt(epsilon) + ux.*sqrt(1+epsilon*ux.^2))/2;
c  = @(ux) sqrt(1+epsilon*ux.^2);
cp = @(ux) epsilon*ux./sqrt(1+epsilon*ux.^2);
fluxes = {{@(ux) -Q(ux), @(ux) -c(ux), @(ux) -cp(ux)}, {Q,c,cp}};

% Initial conditions for L and R. These are output from
%   [L0,R0] = nhyperbolic_characteristic_ICs(A,n,u0,u0t)
% after differentiating
% R0p  = @(x) -2*x.*exp(-x.^2) + 2*epsilon*x.*(3*x.^2 + 1).*exp(-3*x.^2)/9;
% R0pp = @(x) 2*(2*x.^2-1).*exp(-x.^2) + -2*epsilon*x.^2.*(2*x.^2 - 1).*exp(-3*x.^2);
% L0p  = @(x) -2*epsilon*x.^3.*exp(-3*x.^2)/3;
% L0pp = @(x) 2*epsilon*x.^2.*(2*x.^2 - 1).*exp(-3*x.^2);
R0p  = @(x) -x.*exp(-x.^2);
R0pp = @(x) (2*x.^2-1).*exp(-x.^2);
L0p  = R0p;
L0pp = R0pp;
ICs  = {{L0p,L0pp},{R0p,R0pp}};

Nmesh = 100; a = -8; b = 8;
x0 = linspace(a,b,Nmesh); t0 = linspace(0,20,Nmesh);
tic
[X,t,LR,S,tbreaks,partition,XI] = transport_solver_v2(fluxes,ICs,x0,t0,false);
% [X,t,LR,S,tbreaks] = transport_solver(fluxes{1}{1},fluxes{1}{2},fluxes{1}{3},ICs{1}{1},ICs{1}{2},x0,t0,false);
% [X,t,LR,S,tbreaks] = transport_solver(fluxes{2}{1},fluxes{2}{2},fluxes{2}{3},ICs{2}{1},ICs{2}{2},x0,t0,false);
toc

%%
figure
plot(X(1,:),LR(1,:,1),X(end,:),LR(end,:,1))
figure
plot(X(1,:),LR(1,:,2),X(end,:),LR(end,:,2))

% L and R were solved over the same mesh so pointwise operations between
% them are easy to accomplish (wow so expressive!)
% U = LR(1)+LR(2);


