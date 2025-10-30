%% Observe the exact solution of $u_{tt}=A'(u_x)u_{xx}$
% when the two waves do not share common support
clearvars; clc; close all;
config_figures;
this_dir = "interaction_free_dynamics/";
if ~exist(this_dir,'dir')
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

