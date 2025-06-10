clearvars; close all; clc;
% Parameters
L = 8; % x goes from 0 to L (not including the last point)
T = 2; % final time
epsilon = 0.5;
power_of_2 = 4;
base = 2;
h = 0.05;
tau = h/base^power_of_2;

% Discretization
x = 0:h:L;
t = 0:tau:T;
N = length(x);
M = length(t);

% Initial condition: Gaussian at the domain center
u0 = exp(-((x - L/2).^2) / 0.1); % Gaussian
u1 = u0; % Zero initial velocity

% u = CL_conserving_implicit_stencil(epsilon,u0,u1,M,N,tau,h);
u = CL_conserving_semi_implicit_stencil(epsilon,u0,u1,M,N,tau,h);


%% Visualization
x_plot = x(1:5:end);
t_plot = t(1:10:end);
u_plot = u(1:5:end, 1:10:end);

[X_plot, T_plot] = meshgrid(x_plot, t_plot);
figure;
surf(X_plot, T_plot, u_plot');
shading interp; % To smooth the surface
xlabel('x');
ylabel('t');
zlabel('u');
title('Solution of the nonlinear wave equation');

%% Figure 2: Plot of several time snapshots
close all;
figure(2);
hold on;
snapshots = [1, floor(M/4), floor(M/2), floor(3*M/4), M]; % Select several time snapshots
colors = lines(length(snapshots)); % Line colors for different snapshots

for k = 1:length(snapshots)
    plot(x, u(:, snapshots(k)), 'DisplayName', ['t = ' num2str(t(snapshots(k)))], 'Color', colors(k, :), 'LineWidth',1.5);
end
%xlim([5,25]);

xlabel('x');
ylabel('u');
title('Snapshots of u as a function of x');
legend show;
hold off;

%% CLs

clc;

Ut=(u(:,2:end)-u(:,1:end-1))/tau;

Density1 = Ut;
C1 = sum(h*Density1, 1);

Ux =    (u(2:end,1:end-1)-u(1:end-1,1:end-1))/h;
Uxhat = (u(2:end,2:end)-u(1:end-1,2:end))/h;

Density3 = 1/2*(Ux.*Uxhat + Ut(2:end,:).^2)+ 1/12*epsilon*Ux.^2.*Uxhat.^2;
C3 = sum(h*Density3, 1);


Density4 = t(2:end).* Ut(1:end,:) - u(:,2:end);
C4 = sum(h*Density4, 1);



%C4 = C4 - C4(1);
% plot CLs
plotTime = M-1;

figure(3)
clf;


subplot(3,1,1)
plot(t(1:plotTime), C1(1:plotTime), 'linewidth', 1.75);
title('Conservation Laws, 9-pt scheme', 'FontSize', 16, 'interpreter', 'latex');
legend('Density 1', 'fontsize', 15, 'interpreter', 'latex');
hold on;

subplot(3,1,2)
plot(t(1:plotTime), C3(1:plotTime), 'color', [0.9290 0.6940 0.1250], 'linewidth', 1.75);
legend('Density 3', 'fontsize', 15, 'interpreter', 'latex');
%xlabel('Time', 'FontSize', 18, 'interpreter', 'latex');
hold on;

subplot(3,1,3)
plot(t(1:plotTime), C4(1:plotTime), 'color', [1 0 1], 'linewidth', 1.75);
legend('Density 4', 'fontsize', 15, 'interpreter', 'latex');
xlabel('Time', 'FontSize', 18, 'interpreter', 'latex');
hold off;

%% CL 2 not conserved


figure(22)
clf;

Density2 = Ux.*Ut(2:end,:);
C2 = sum(h*Density2, 1);
plot(t(1:plotTime), C2(1:plotTime), 'color', [1 0 1], 'linewidth', 1.75);
legend('Density 2', 'fontsize', 15, 'interpreter', 'latex');
xlabel('$t$', 'FontSize', 18, 'interpreter', 'latex');

