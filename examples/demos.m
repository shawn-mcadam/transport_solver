%% Test script for Burgers' solver
clearvars; clc; close all;
config_figures
this_dir = "demo/";
breakcolor = "#A2142F";
coloru  = "#53A966";
colorux  = "#c3e166";

if exist(this_dir,'dir') ~= 7
    mkdir(this_dir)
end


%% --- Solve a common test for Burger's equation: u_t + u*u_x = 0 ---
% The flux
Q = @(u) u.^2/2; c = @(u) u; cp = @(u) 1;
% Initial condition and its derivative
phi = @(x) exp(-x.^2); phip = @(x) -2*x.*exp(-x.^2);

% fine mesh of single valued solution
Nmesh = 350; a = -6; b = 8; tf = 14;
x0 = linspace(a,b,Nmesh); t0 = linspace(0,tf,Nmesh);
[X,t,U,S,tbreaks] = transport_solver(Q,c,cp,phi,phip,x0,t0,false);
idxs = ~isnan(X(1,:)); Ijump = find(isnan(X(1,:)));

%% compute the characteristics
chars = XI(end,idxs) + c(r0(XI(end,idxs))).*t;
for i=2:15:length(t)
    chars(1:i,end+1:end+2*numel(Ijump)) = XI(i,[Ijump-1,Ijump+1]) + c(r0(XI(i,[Ijump-1,Ijump+1]))).*t(1:i);
    chars(i+1:end,end-2*numel(Ijump):end) = NaN;
end
chars = sortrows(chars')';


f = figure; f.Position = [680 458 900 260];
pause(1e-5)
mindelay = 0.025;
gifname = 'burgers_gaussian.gif';
for k = 1:length(t)
    tiledlayout
    nexttile
    pcolor(X,t,U,EdgeColor='none'); hold on
    % plot(X,t,color="black",LineWidth=0.5); % computational mesh
    plot(chars,t,color="black",LineWidth=0.5); % characteristics mesh
    plot(S,t,"-",color=breakcolor,LineWidth=2);
    yline(tbreaks);
    colormap summer
    colorbar
    xlabel("$x$"); ylabel("$t$")
    xlim tight; ylim tight
    hold off

    nexttile
    plot(X(k,:), U(k,:),color=coloru); hold on
    for i = find(tbreaks < t(k))
    plot([X(k,Ijump(i)-1);X(k,Ijump(i)+1)], [U(k,Ijump(i)-1);U(k,Ijump(i)+1)],'--',color=[0.7,0.7,0.7])
    end
    hold off

    xlim tight; ylim([0,1])
    title("$t="+string(round(t(k),2))+"$")
    legend("$u$")

    if k == 1, gif(gifname,'nodither','DelayTime',mindelay); else, gif; end
end
movefile(gifname, this_dir);



%% Repeat the demo but with the nonlinear string equation
syms Sx Sepsilon Sux

% Moving frame of reference removes the O(1) translation over time
A(Sux) = atanh(Sux/sqrt(1+Sux^2));
Sc(Sux) = simplify(sqrt(diff(A(Sux),Sux)))-1;

% Derive the initial conditions
R0(Sx)  = exp(-0.3*Sx^2);
L0(Sx)  = 0*Sx;

% define flux functions Q,c,cp
Q  = matlabFunction(int(Sc,Sux));
c  = matlabFunction(Sc);
cp = matlabFunction(diff(Sc,Sux));
r0 = matlabFunction(diff(R0,Sx));
r0p = matlabFunction(diff(R0,Sx,Sx));


% fine mesh of single valued solution
Nmesh = 100; a = -7; b = 7; tf = 150;
x0 = linspace(a,b,Nmesh); t0 = linspace(0,tf,3*ceil(tf));
[X,t,Ux,S,tbreaks,~,XI] = transport_solver(Q,c,cp,r0,r0p,x0,t0,false);
idxs = ~isnan(X(1,:)); Ijump = find(isnan(X(1,:)));
U = NaN*ones(size(Ux));
for i = 1:length(t)
    U(i,idxs) = cumtrapz(X(i,idxs),Ux(i,idxs));
end

%% Compute the characteristics
chars = XI(end,idxs) + c(r0(XI(end,idxs))).*t;
for i=2:15:length(t)
    chars(1:i,end+1:end+2*numel(Ijump)) = XI(i,[Ijump-1,Ijump+1]) + c(r0(XI(i,[Ijump-1,Ijump+1]))).*t(1:i);
    chars(i+1:end,end-2*numel(Ijump):end) = NaN;
end
chars = sortrows(chars')';

%%
pcolor(X,t,U,EdgeColor='none'); hold on
plot(chars,t,color="black",LineWidth=0.5);
plot(S,t,"-",color=breakcolor,LineWidth=2);
yline(tbreaks);
colormap summer
colorbar
xlabel("$\tilde x = x - t$"); ylabel("$t$")
legend("$u$")
xlim tight; ylim tight
hold off

%% f = figure; f.Position = [680 458 1100 350];
f = figure; f.Position = [680 458 1200 500];
pause(1e-5)
mindelay = 0.025;
% gifname = 'string_gaussian.gif';

v = VideoWriter("string_gaussian.avi",'MPEG-4');
% v.Quality=100;
open(v);
% for k = 1:length(t)
for k = 1:50
    plot(X(k,:)+t(k), U(k,:),color=coloru); hold on
    plot(X(k,:)+t(k), Ux(k,:),color=colorux);
    for i = find(tbreaks < t(k))
    plot([X(k,Ijump(i)-1);X(k,Ijump(i)+1)]+t(k), [Ux(k,Ijump(i)-1);Ux(k,Ijump(i)+1)],'--',color=[0.7,0.7,0.7])
    end
    hold off
    
    xlim tight; ylim([-0.5,1]);
    xlabel("$x$")
    title("$t="+string(round(t(k),2))+"$")
    legend("$u$","$u_x$")

    frame = getframe(gcf);
    writeVideo(v,frame)
end
close(v);