%Radiation/chemical leak (running example) model parameter estimation
clc; clear; close all;
%% Adjustable Parameters
% General simulation parameters
M = 3;                  % Number of robots
dt = 0.5;               % Sampling period [s]
sim_time = 15;          % Simulation time [s]
K = sim_time/dt;        % Total # of simulation steps
Ts = 0.5;                 % MPC sampling period
sim_params = [Ts, dt];
state_plot_selec = [1];            % Select states to plot
error_cv_selec = [1,3,5,7];
sc = [200 1 0.05 0.0001 20 0.1 20 0.1]';        % Scaling (FIXED)
%sc = [200 1 0.05 1 20 1 20 1]';
%sc = [1 1 1 1 1 1 1 1]';
T = diag(sc);

% Random seed
rng(57)

% MPC parameters
Hp = 3;             % Prediction horizon
Hu = 3;             % Control horizon
mpc_params = [Hp, Hu];           % Hp, Hu
cost_params = [0.1, 1, 1, 1];    % lambda1, lambda2, epsilon
min_dist = 1;                  % Minimum distance between robots

% Map bounderies
xmin = 0; xmax = 40;
ymin = 0; ymax = 40;
zmin = 0; zmax = 5;
map_bounds = [xmin, xmax, ymin, ymax];

% True initial process states
M_ref = 500;
beta_ref = 0.05;
xs_ref = 20;
ys_ref = 20;
M_dot_ref = 0;
beta_dot_ref = 0;
xs_dot_ref = 0;
ys_dot_ref = 0;
z_0 = [M_ref; M_dot_ref; beta_ref; beta_dot_ref; xs_ref; xs_dot_ref; ys_ref; ys_dot_ref]; % Initial true state vector
Nx_p = length(z_0); % Number of states in process state vector

% Guessed initial process states
M_0 = 400;
beta_0 = 0.05;
xs_0 = 15;
ys_0 = 25;
M_dot_0 = 0;
beta_dot_0 = 0;
xs_dot_0 = 0;
ys_dot_0 = 0;
z_est_0 = [M_0; M_dot_0; beta_0; beta_dot_0; xs_0; xs_dot_0; ys_0; ys_dot_0]; % Initial guessed state vector

% Initial error covariance matrix
P_0 = zeros(Nx_p);
P_0(1,1) = (M_ref - M_0)^2;
P_0(3,3) = (beta_ref - beta_0)^2;
P_0(5,5) = (xs_ref - xs_0)^2;
P_0(7,7) = (ys_ref - ys_0)^2;
%P_0 = eye(Nx_p);

% Initial robot positions, Z0 and control inputs, Uprev
x_1 = zeros(2*M,1); % (assuming M=4, so 2*M values)
for m = 1:M
    x_1(2*m-1) = m;
    x_1(2*m) = 0;
end
U_prev = zeros(2*M,1);           % (2*M values again)

% MPC Tuning
Q_vec = [1/z_est_0(1); 0; 1/z_est_0(3); 0; 1/z_est_0(5); 0; 1/z_est_0(7); 0];
R_single = eye(2);
mpc_tuning = {Q_vec, R_single};

%% Dynamic Parameter Model
tau_beta = 0.97;       %spread increase parameter (exponential decay of beta)

% Inputs
Nu_p = 2;
beta0 = 0.001; %steady-state value for beta
u_beta = beta0*(1-tau_beta)/dt; % constant input to achieve beta0 in steady-state
u_M = 5; % Leakage rate offset
u = ones(Nu_p,K+Hp+1);
u(1,:) = u(1,:)*u_M;
u(2,:) = u(2,:)*u_beta;

% Process system matrix (depends on dt)
A_func = @(dt) [1 dt 0 0 0 0 0 0;
           0 1 0 0 0 0 0 0;
           0 0 tau_beta dt 0 0 0 0;
           0 0 0 1 0 0 0 0;
           0 0 0 0 1 dt 0 0;
           0 0 0 0 0 1 0 0;
           0 0 0 0 0 0 1 dt;
           0 0 0 0 0 0 0 1];
B_func = @(dt) [dt 0;
           0 0;
           0 dt;
           0 0;
           0 0;
           0 0;
           0 0;
           0 0];

% Create system matrices with corresponding sampling periods
A = A_func(dt);
B = B_func(dt);
A_sc = inv(T)*A*T;  % Scaled
B_sc = inv(T)*B;    % Scaled

% Process input matrix
% B = zeros(Nx_p,2);
% B(1,1) = dt; B(3,2) = dt;

% Process noise (Q)
G_func = @(dt) [0.5*dt^2 0 0 0;
    dt 0 0 0;
    0 0.5*dt^2 0 0;
    0 dt 0 0;
    0 0 0.5*dt^2 0;
    0 0 dt 0;
    0 0 0 0.5*dt^2;
    0 0 0 dt];
% 
G = G_func(dt);
sigma_M = 2; sigma_beta = 0.00001; sigma_x = 0.05; sigma_y = 0.05;
v = [sigma_M sigma_beta sigma_x sigma_y]';
mu_w = zeros(length(v),1);
Q = (inv(T)*G)*diag([sigma_M^2 sigma_beta^2 sigma_x^2 sigma_y^2])*(inv(T)*G)'; % This is the process noise covariance matrix

% Measurement noise
sigma_measurement = 0.1;
mu_r = zeros(M,1);
R = sigma_measurement^2*eye(M); % measurement noise covariance

% Everyting recomputed for the MPC
A_MPC = A_func(Ts);
B_MPC = B_func(Ts);
Q_MPC = G_func(Ts)*diag([sigma_M^2 sigma_beta^2 sigma_x^2 sigma_y^2])*G_func(Ts)';
Q_MPC_sc = (inv(T)*G_func(Ts))*diag([sigma_M^2 sigma_beta^2 sigma_x^2 sigma_y^2])*(inv(T)*G_func(Ts))';
A_MPC_sc = inv(T)*A_MPC*T;  % Scaled
B_MPC_sc = inv(T)*B_MPC;    % Scaled

%% Robot State Space Model
% Integrator
A_r_single_func = @(dt) [1 0;       % Function to compute robot system matrix from dt
                         0 1];
B_r_single_func = @(dt) [dt 0;      % Function to compute robot input matrix from dt
                         0 dt];

% Compute robot A_r and B_r (simulation and MPC are different)
A_r_single = A_r_single_func(dt); % z = [x;y]
B_r_single = B_r_single_func(dt); % u = [vx;vy]
A_r_single_MPC = A_r_single_func(Ts);
B_r_single_MPC = B_r_single_func(Ts);

% Construct full robot dynamics matrix
A_r = kron(eye(M), A_r_single);
B_r = kron(eye(M), B_r_single);
A_r_MPC = kron(eye(M), A_r_single_MPC);
B_r_MPC = kron(eye(M), B_r_single_MPC);

% Save number of states and inputs for robot model
Nx_r = size(A_r_single,1);
Nu_r = size(B_r_single,2);

%% Energy
% Energy variables
e = zeros(M,K+1);
e(:,1) = ones(M,1)*0.8;     % Robots fully charged at time 0
charge_control = ones(M,K+1);
charge_control(:,1) = zeros(M,1);

% Charging params
h_ch = 5; % Size of station
s_ch = 0.5; % Sharpness of corners
cx = 0; % x-coordinate of station center
cy = 0; % y-coordinate of station center

% Energy dynamics
charge_rate = 0.2;
%en_charge = @(x,y) (1/(1+exp(t1*(x-o1))))*(1/(1+exp(t1*(y-o1))))*charge_rate; % Energy charging function
en_charge = @(x,y) (1/(1+exp(-s_ch*(x-(cx-h_ch)))))*(1/(1+exp(s_ch*(x-(cx+h_ch)))))*(1/(1+exp(-s_ch*(y-(cy-h_ch)))))*(1/(1+exp(s_ch*(y-(cy+h_ch)))))*charge_rate; 
en_cons = @(v) 0.0001*v^2 + 0.05; % Energy consumption function

%% Create Required Vectors and Matrices
% MPC
x = zeros(2*M,K+1);         % Robot positions

% KALMAN FILTER           
z = zeros(Nx_p,K+1);           % True process states (Nx_p)
z_est = zeros(Nx_p,K+1);       % Estimated process states by Kalman filter
P = zeros(Nx_p,Nx_p,K+1);         % Error covariance matrix
y = zeros(M,K+1);           % Measurements of process

% Put initial conditions into vectors
z(:,1) = inv(T)*z_0;
z_est(:,1) = inv(T)*z_est_0;
x(:,1) = x_1;
x(:,2) = x_1;
P(:,:,1) = P_0;

%% Plotting preparation
% Initializations
x_axis = linspace(xmin,xmax,50);    % Process x-axis
y_axis = linspace(ymin,ymax,50)';   % Process y-axis
z_values_true = get_h(z(:,1),x_axis,y_axis,sc);    % Process values in defined area
z_values_est = get_h(z_est(:,1), x_axis, y_axis,sc);
t_vec = (0:K) * dt;               % Time vector for plotting
colors = {'b','r','g','k','y','c','m','b'};     % valid MATLAB colors

% Initialize plotting objects that can be updated 
process_plot_true = gobjects(1,1);         
process_plot_est = gobjects(1,1); 

% Initial process plot setup
figure(1)
set(gcf, 'Position',  [100, 300, 1400, 400])
%  True process plot
subplot(1,3,1)
hold on;
axis([xmin xmax ymin ymax zmin zmax]);
xlabel('X'); ylabel('Y');
process_plot_true(1) = mesh(x_axis,y_axis,z_values_true, 'FaceAlpha','0.7', 'EdgeAlpha','0.7');
view([45 45])
% Estimated process plot
subplot(1,3,2)
hold on;
axis([xmin xmax ymin ymax zmin zmax]);
xlabel('X'); ylabel('Y');
process_plot_est(1) = mesh(x_axis,y_axis,z_values_est, 'FaceAlpha','0.7', 'EdgeAlpha','0.7');
view([45 45])

% Store plot handles
ax1 = subplot(1,3,1); hold(ax1,'on');
ax2 = subplot(1,3,2); hold(ax2,'on');

% Create robot quivers (they are identical, but must be created twice to live in multiple plots)
robot_true = gobjects(M,1);     
robot_est = gobjects(M,1);

% Plot robots on process plot (with their inital conditions)
subplot(1,3,1)
thickness = 5;
for m = 1:M
    robot_true(m) = plot(ax1, x_1(2*m-1), x_1(2*m), 'o', ...
                                    'MarkerSize', thickness, ...
                                    'MarkerFaceColor', colors{m}, ...
                                    'MarkerEdgeColor', 'k', ...
                                    'LineWidth', 1);
    robot_est(m) = plot(ax2, x_1(2*m-1), x_1(2*m), 'o', ...
                                    'MarkerSize', thickness, ...
                                    'MarkerFaceColor', colors{m}, ...
                                    'MarkerEdgeColor', 'k', ...
                                    'LineWidth', 1);
end



% Plot selected states
z_selected = zeros(length(state_plot_selec),K+1);          % Vector for storage of selected true process states
z_est_selected = zeros(length(state_plot_selec),K+1);      % Vector for storage of selected estimated process states
z_selected_plot = gobjects(size(z_selected,1),1);
z_est_selected_plot = gobjects(size(z_selected,1),1);
subplot(1,3,3)
hold on
for i = 1:length(state_plot_selec)
    s = state_plot_selec(i);
    z_selected(i,1) = T(s,s)*z(s,1);   % Plug in initial states as first column
    z_est_selected(i,1) = T(s,s)*z_est(s,1);   % Plug in initial guessed states as first column
    z_selected_plot(i) = plot(0,z_selected(1,1),'DisplayName',sprintf('State %d', s),'Color',colors{i});
    z_est_selected_plot(i) = plot(0,z_est_selected(1,1),'DisplayName',sprintf('State %d est.', s),'Color',colors{i},'LineStyle','--');
end

% Plot a selection of true and guessed initial states of the process
xlabel('Time [s]')
legend('Interpreter','latex','Location','west')

%% Simulation loop
x_opt = zeros(2*M,K+1);         % Store optimal robot positions
u_opt = U_prev;
% rob_init = [x_opt, u_opt];
process_params = zeros(4, K+1);

% for random movement
robot_theta = rand(M,1) * 2*pi;
robot_v = zeros(M,1);          % Linear velocities
robot_omega = zeros(M,1);      % Angular velocities
max_speed = 0.5*2;        % Maximum robot speed
max_omega = pi/2;       % Maximum angular velocity
change_prob = 0.1;      % Probability of changing direction

z_hat = zeros(Nx_p,K+1);
z_hat(:,1) = z_est(:,1);

% Warm start
sol_prev = 0;
do_warm_start = 0;

% Live plotting variables
z_real = zeros(Nx_p,K+1);
z_est_real = zeros(Nx_p,K+1);
z_real(:,1) = z_0;
z_est_real(:,1) = z_est_0;

% Offline plotting varialbes
e_low_slack = zeros(M,K+1);     % Slack on lower energy bound
e_high_slack = zeros(M,K+1);    % Slack on higher energy bound
b_slack = zeros(M,K+1);         % Slack on energy barrier
dist_slack = zeros(M,K+1);      % Slack on minimum collision distance

% MAIN SIMULATION LOOP
for k=2:K+1
    % Update true process
    w = (inv(T)*G)*normrnd(mu_w,v);
    z(:,k) = A_sc*z(:,k-1) + B_sc*u(:,k-1) + w;
    
    % Update robot positions
    x(:,k) = A_r*x(:,k-1) + B_r*u_opt;

    % Energy dynamics
    for m = 1:M
        % Check if the robot were at the charging station at the time step before
        % if((x(2*m-1,k-1) - charger_x)^2 + (x(2*m,k-1) - charger_y)^2 <= charger_r^2)
        %     charging(m,k-1) = 1;
        % else
        %     charging(m,k-1) = 0;
        % end
        % if(e(m,k-1) >= 1-charge_rate*dt)
        %     e_full(m,k-1) = 1;
        % else
        %     e_full(m,k-1) = 0;
        % end
        % Compute new energy
        e(m,k) = e(m,k-1) + (charge_control(m,k-1)*en_charge(x(2*m-1,k), x(2*m,k)) - en_cons(sqrt(u_opt(2*m-1)^2 + u_opt(2*m)^2)))*dt;
    end
    %disp(e(:,k));

    % Take measurements
    r = mvnrnd(mu_r',R)';
    y(:,k) = get_h_vec(z(:,k), x(:,k), M, sc) + r;
    
    % Run Kalman filter iteration
    if(k ==180)
        disp("")
    end
    [z_est(:,k), P(:,:,k), z_hat(:,k)] = EKF(z_est(:,k-1), P(:,:,k-1), A_sc, B_sc, u(:,k-1), Q, y(:,k), R, x(:,k), sc);
   
    % Package all Kalman filter information for MPC
    u_mpc = u(:,k:k+Hp);                % This could be avoided by removing DeltaT in B and in the u_beta calc.
    u_mpc(2,:) = u_mpc(2,:)*(dt/Ts);    % Scale input to match MPC sampling period
    KF_params = {z_est(:,k), P(:,:,k), A_MPC_sc, B_MPC_sc, u_mpc, Q_MPC_sc, R, error_cv_selec};

    % Package all robot information
    ROB_params = {x(:,k), A_r_MPC, B_r_MPC, u_opt};

    % Package energy information
    EN_params = {e(:,k), en_charge, en_cons, cx, cy};
    
    % Compute optimal stuff
    [X_opt, U_opt, P_trace, sol_prev, charge_control(:,k), slacks] = MPC_func(ROB_params, KF_params, EN_params, mpc_tuning, mpc_params, cost_params, M, map_bounds, min_dist, sim_params, 2, do_warm_start, sol_prev, sc);
    u_opt = U_opt(:,2);

    % set warm start to 1 after first iteration
    do_warm_start = 1; 
    
    % Insert slack variables into respective variables for later plotting
    e_low_slack(:,k+1) = slacks{1};
    e_high_slack(:,k+1) = slacks{2};
    b_slack(:,k+1) = slacks{3};
    dist_slack(:,k+1) = slacks{4};

    % Update process plot
    %subplot(1,2,1)
    z_values_true = get_h(z(:,k),x_axis,y_axis, sc);
    z_values_est = get_h(z_est(:,k),x_axis,y_axis, sc);
    set(process_plot_true(1), 'XData', x_axis, 'YData', y_axis, 'ZData', z_values_true)
    set(process_plot_est(1), 'XData', x_axis, 'YData', y_axis, 'ZData', z_values_est)

    % Update robot positions
    for m = 1:M
        set(robot_true(m),'XData', x(2*m-1,k),'YData', x(2*m,k));
        set(robot_est(m),'XData', x(2*m-1,k),'YData', x(2*m,k));
    end
    
    %update state and estimate plot
    z_real(:,k) = T*z(:,k);
    z_est_real(:,k) = T*z_est(:,k);
    for i = 1:length(state_plot_selec)
        s = state_plot_selec(i);
        z_selected(i,k) = z_real(s,k);
        z_est_selected(i,k) = z_est_real(s,k);
        set(z_selected_plot(i), 'XData', t_vec(1:k),'YData', z_selected(i,1:k))
        set(z_est_selected_plot(i), 'XData', t_vec(1:k), 'YData', z_est_selected(i,1:k))
    end
    drawnow limitrate
end

% Cut off last column of slack variables before plotting
e_low_slack(:,end) = [];
e_high_slack(:,end) = [];
b_slack(:,end) = [];
dist_slack(:,end) = [];

%% Plotting After Finish
close all

% Figure settings
hFig = figure; 

set(hFig, ...
    'Units','centimeters', ...
    'Position',[5 5 20 14] ...     % [left bottom width height] on screen
);

set(hFig, ...
    'PaperSize',[21 29.7], ...
    'PaperOrientation','portrait', ...
    'PaperPosition',[1 1 19 27.7], ...
    'PaperPositionMode','manual' ...
);

set(hFig, ...
    'DefaultAxesFontSize', 10, ...
    'DefaultTextFontSize', 10, ...
    'DefaultLegendFontSize', 9 ...
);

% Time vector
t = t_vec;
subplot(4,2,1); hold on
hold on
colors = lines(M);
for m = 1:M
    plot(x(2*m-1, :), x(2*m, :), '-o', 'MarkerSize', 1.5)
    plot(x(2*m-1, 1), x(2*m, 1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 3, 'HandleVisibility', 'off')
    plot(x(2*m-1, end), x(2*m, end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3, 'HandleVisibility', 'off')
end
xlabel('X [m]')
ylabel('Y [m]')
title('Robot Trajectories (Top View)')
axis([xmin xmax ymin ymax])
% legend('Location','best')
grid on
hold off

% Energy
subplot(4,2,2)
hold on 
for m = 1:M
    plot(t, e(m,:), 'LineWidth',0.9)
end
xlabel('Time [s]')
ylabel(sprintf('Percentage'))
title(sprintf('Robot Energies'))
ylim([0 1])
grid on
hold off

% Error covariance
subplot(4,2,7)
hold on 
semilogy(t, squeeze(P(1,1,:)), 'LineWidth',0.9)
semilogy(t, squeeze(P(3,3,:))+1e-3, 'LineWidth',0.9)
semilogy(t, squeeze(P(5,5,:)), 'LineWidth',0.9)
semilogy(t, squeeze(P(7,7,:)), 'LineWidth',0.9)
xlabel('Time [s]')
ylabel(sprintf('Value'))
title(sprintf('State Variances'))
grid on
hold off
set(gca,'YScale','log')

% Slacks
subplot(4,2,8)
hold on 
for m = 1:M
    plot(t, e_low_slack(m,:), 'LineWidth',0.9)
    plot(t, e_high_slack(m,:), 'LineWidth',0.9)
    plot(t, b_slack(m,:), 'LineWidth',0.9)
    plot(t, dist_slack(m,:), 'LineWidth',0.9)
end
xlabel('Time [s]')
ylabel(sprintf('Value'))
title(sprintf('Slack Variables'))
grid on
hold off

% State plot of M
subplot(4,2,3)
plot( t, z_real(1, :),   'b-', 'LineWidth',0.9, 'DisplayName','True' )
hold on
plot( t, z_est_real(1,:), 'r--','LineWidth',0.9, 'DisplayName','Estimated' )
xlabel('Time [s]')
ylabel(sprintf('M'))
title(sprintf('Process State M'))
%legend('Location','best')
ylim([0 1000])
grid on
hold off

% State plot of beta
subplot(4,2,4)
plot( t, z_real(3, :),   'b-', 'LineWidth',0.9, 'DisplayName','True' )
hold on
plot( t, z_est_real(3,:), 'r--','LineWidth',0.9, 'DisplayName','Estimated' )
xlabel('Time [s]')
ylabel(sprintf('\\beta'))
title(sprintf('Process State \\beta'))
% legend('Location','best')
grid on
hold off

% State plot of xs
subplot(4,2,5)
plot( t, z_real(5, :),   'b-', 'LineWidth',0.9, 'DisplayName','True' )
hold on
plot( t, z_est_real(5,:), 'r--','LineWidth',0.9, 'DisplayName','Estimated' )
xlabel('Time [s]')
ylabel(sprintf('x_s'))
title(sprintf('Process State x_s'))
% legend('Location','best')
ylim([xmin xmax])
grid on
hold off

% State plot of ys
subplot(4,2,6)
plot( t, z_real(7, :),   'b-', 'LineWidth',0.9, 'DisplayName','True' )
hold on
plot( t, z_est_real(7,:), 'r--','LineWidth',0.9, 'DisplayName','Estimated' )
xlabel('Time [s]')
ylabel(sprintf('y_s'))
title(sprintf('Process State y_s'))
% legend('Location','best')
ylim([ymin ymax])
grid on
hold off

print(hFig, 'myFigure.pdf', '-dpdf', '-bestfit');

% hFig = figure;
% set(hFig, 'Units','centimeters', 'Position',[5 5 20 14], ...
%           'PaperSize',[21 29.7], 'PaperOrientation','portrait', ...
%           'PaperPosition',[1 1 19 27.7], 'PaperPositionMode','manual', ...
%           'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, ...
%           'DefaultLegendFontSize',9);
% 
% t = tiledlayout(6,4, ...
%      'TileSpacing','compact', ...   % remove gaps between tiles
%      'Padding','none');             % remove outer margins
% 
% ax1 = nexttile(1, [2 2]);
% plot(ax1, t_vec, z_real(1,:),   'b-', 'LineWidth',0.9, 'DisplayName','True');
% hold(ax1,'on');
% plot(ax1, t_vec, z_est_real(1,:),'r--','LineWidth',0.9, 'DisplayName','Estimated');
% hold(ax1,'off');
% title(ax1,'Process State M');
% xlabel(ax1,'Time [s]');
% ylabel(ax1,'M');
% ylim(ax1,[0 1000]);
% legend(ax1,'Location','best');
% grid(ax1,'on');
% 
% ax2 = nexttile(3, [2 2]);
% plot(ax2, t_vec, z_real(3,:),   'b-', 'LineWidth',0.9,'DisplayName','True');
% hold(ax2,'on');
% plot(ax2, t_vec, z_est_real(3,:),'r--','LineWidth',0.9,'DisplayName','Estimated');
% hold(ax2,'off');
% title(ax2,'Process State \beta');
% xlabel(ax2,'Time [s]');
% ylabel(ax2,'\beta');
% legend(ax2,'Location','best');
% grid(ax2,'on');
% 
% ax3 = nexttile(9, [2 2]);
% plot(ax3, t_vec, z_real(5,:),   'b-', 'LineWidth',0.9,'DisplayName','True');
% hold(ax3,'on');
% plot(ax3, t_vec, z_est_real(5,:),'r--','LineWidth',0.9,'DisplayName','Estimated');
% hold(ax3,'off');
% title(ax3,'Process State x_s');
% xlabel(ax3,'Time [s]');
% ylabel(ax3,'x_s');
% ylim(ax3,[xmin xmax]);
% legend(ax3,'Location','best');
% grid(ax3,'on');
% 
% ax4 = nexttile(11, [2 2]);
% plot(ax4, t_vec, z_real(7,:),   'b-', 'LineWidth',0.9,'DisplayName','True');
% hold(ax4,'on');
% plot(ax4, t_vec, z_est_real(7,:),'r--','LineWidth',0.9,'DisplayName','Estimated');
% hold(ax4,'off');
% title(ax4,'Process State y_s');
% xlabel(ax4,'Time [s]');
% ylabel(ax4,'y_s');
% ylim(ax4,[ymin ymax]);
% legend(ax4,'Location','best');
% grid(ax4,'on');
% 
% ax5 = nexttile(18, [2 2]);
% hold(ax5,'on');
% colors = lines(M);
% for m = 1:M
%     plot(ax5, x(2*m-1,:), x(2*m,:), '-o', 'Color',colors(m,:), 'LineWidth',1);
%     plot(ax5, x(2*m-1,1), x(2*m,1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 3, 'HandleVisibility', 'off');
%     plot(ax5, x(2*m-1,end), x(2*m,end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3, 'HandleVisibility', 'off');
% end
% hold(ax5,'off');
% axis(ax5,[xmin xmax ymin ymax]);
% xlabel(ax5,'X [m]');
% ylabel(ax5,'Y [m]');
% title(ax5,'Robot Trajectories (Top View)');
% legend(ax5,'Location','best');
% grid(ax5,'on');
% 
% print(hFig,'myFigure.pdf','-dpdf','-bestfit');
