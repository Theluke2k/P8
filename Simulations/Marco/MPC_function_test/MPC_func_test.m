clc; clear; close all;
% addpath('C:\Users\marco\OneDrive - Aalborg Universitet\Documents\MATLAB\casadi-3.7.0-windows64-matlab2018b')

% MPC and cost params
Hp = 6; Hu = 3;
mpc_params = [Hp, Hu];             % Hp, Hu
cost_params = [1, 1e4, 1e-6];    % lambda1, lambda2, epsilon

% Environment
xmin = 0; xmax = 20;
ymin = 0; ymax = 20;
map_bounds = [xmin, xmax, ymin, ymax];     % xmin, xmax, ymin, ymax
N = 3;
min_dist = 0.5;                    % Minimum distance between robots
sim_params = [1, 0.1];           % Ts, dt

% Example initialization values
Z0 = zeros(2*N,1); % (assuming N=4, so 2*N values)
for j = 1:N
    Z0(2*j-1) = j;
end
Uprev = zeros(2*N,1);           % (2*N values again)
rob_init = [Z0, Uprev];

Z0_KF = [1e-3; 0; 1e-3; 0; 8; 0; 12; 0];
P0_KF = 1e-1*ones(8);
for i = 1:2:3
    P0_KF(i,i) = 1;
end
kf_init = [Z0_KF, P0_KF];

% Call the function
[Z_sol, U_sol, P_trace] = MPC_object(rob_init, kf_init, mpc_params, cost_params, N, map_bounds, min_dist, sim_params, 3);

%% PLOTS
figure;
subplot(1,2,1);
plot(1:Hp+1, P_trace, '-o');
xlabel("Timestep (k)");
ylabel("tr(P)");
grid on;

subplot(1,2,2);
hold on;
for j = 1:N
    plot(Z_sol(2*j-1,:), Z_sol(2*j,:), '-o', DisplayName=sprintf("Robot %d",j));
    plot(Z_sol(2*j-1,1), Z_sol(2*j,1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 3, 'HandleVisibility', 'off');
    plot(Z_sol(2*j-1,end), Z_sol(2*j,end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3, 'HandleVisibility', 'off');
end
hold off;
legend;
xlabel("x");
ylabel("y");
xlim([xmin, xmax]);
ylim([ymin ymax]);
grid on;


