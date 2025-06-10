clear all;
close all;
clc;

% Parameters
% L = 10; % Length of the spatial domain
% T = 5;  % Total time
% h = 0.1; % Spatial step size
% tau = 0.01; % Time step size

L = 3; % x goes from -L to L (not inclusive of last point)
T = 2; % final time

epsilon = 0.2;

power_of_2 = 4;%2
base = 2;

h=0.035;%0.05 fast; 0.01 doable but long
tau = h/base^power_of_2;


% Discretization
x = 0:h:L;
t = 0:tau:T;
N = length(x);
M = length(t);

% Initial condition: Gaussian at the domain center
u0 = exp(-((x - L/2).^2) / 0.1); % Gaussian
u1 = u0  ; % Zero initial velocity

% Initialize the solution matrix
u = zeros(N, M);
u(:, 1) = u0;
u(:, 2) = u1;

% Periodic boundary conditions
u(1, :) = u(end-1, :); % because x(1) == x(N)
u(end, :) = u(2, :);

% % Solver for the implicit scheme using fzero
% for n = 2:M-1
%     u_prev = u(:, n-1);
%     u_curr = u(:, n);
% 
%     for i = 2:N-1
%         %func = @(u_next_i) u_next_i - tau^2 * (1 + ((u_curr(i+1) - u_curr(i-1)) / (2*h))^2) * (u_curr(i+1) - 2*u_next_i + u_curr(i-1)) / h^2 - 2*u_curr(i) + u_prev(i);
%         func = @(u_next_i) u_next_i - tau^2 * epsilon * (1 + ((u_curr(i+1) - u_curr(i-1)) / (2*h))^2) * (u_curr(i+1) - 2*u_curr(i) + u_curr(i-1)) / h^2 - 2*u_curr(i) + u_prev(i);
%         u(i, n+1) = fzero(func, u_curr(i));
%     end
% 
%     % Enforce periodic boundary conditions
%     u(1, n+1) = u(N-1, n+1);
%     u(N, n+1) = u(2, n+1);
% end

% Solver for the implicit scheme using fsolve
for n = 2:M-1
    u_prev = u(:, n-1);
    u_curr = u(:, n);
    
    % calculate only for indices 2..end-1. Then 1<-end-1; end<-2
    %still do comp of N elemens, then impose periodicity

    uxx = (u_curr([2:end 1]) - 2*u_curr + u_curr([end 1:end-1])) / h^2;

    ux2 =  ( u_curr([2:end 1]) - u_curr ).^2/h^2; 
    ux2m = ( u_curr - u_curr([end 1:end-1])) .^2/h^2; 

    Unextx = @(unext) ( unext([2:end 1]) - unext )/h;
    Unextxm = @(unext) ( unext - unext([end 1:end-1]) )/h;

    Uprevx = ( u_prev([2:end 1]) - u_prev )/h; 
    Uprevxm = ( u_prev - u_prev([end 1:end-1]) )/h; 

%    F = @(u_next) u_next- 2*u_curr + u_prev - tau^2 * uxx -1/6*tau^2.*...
%          ( ux2.*(Unextx(u_next)+Uprevx)-ux2m.*(Unextxm(u_next)+Uprevxm) )/h;

   F = @(u_next) (u_next- 2*u_curr + u_prev)/tau^2 - uxx -1/6*...
        epsilon* ( ux2.*( (( u_next([2:end 1]) - u_next )/h) +Uprevx)-ux2m.*( (( u_next - u_next([end 1:end-1]) )/h)  +Uprevxm) )/h;

  
    % Solve the nonlinear system using fsolve
    options = optimoptions('fsolve', 'Display', 'none', 'TolFun', 1e-6, 'TolX', 1e-6);
    u_next = fsolve(F, u_curr, options);
    
    % Enforcing periodic boundary conditions for u_next
    u(:, n+1) = u_next;
    % u(1, n+1) = u(end-1, n+1);
    % u(end, n+1) = u(2, n+1);
end

%% Visualization
close all;
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

%u(space,time)
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


%figure(3);
%plot(t(1:end-1),C3, 'linewidth', 1.75);

%C4 = C4 - C4(1);
%% plot CLs
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

