%Radiation/chemical leak (running example) model parameter estimation
clc; clear; close all;
%% Adjustable Parameters
% General simulation parameters
M = 3;                  % Number of robots
dt = 0.5;               % Sampling period [s]
sim_time = 100;          % Simulation time [s]
K = sim_time/dt;        % Total # of simulation steps
Ts = 0.5;                 % MPC sampling period
sim_params = [Ts, dt];
state_plot_selec = [1];            % Select states to plot
error_cv_selec = [1,3,5,7];
%sc = [200 1 0.05 0.0001 20 0.1 20 0.1]';        % Scaling
sc = [1 1 1 1 1 1 1 1]';
T = diag(sc);

% Random seed
rng(51)

% MPC parameters
Hp = 6;             % Prediction horizon
Hu = 6;             % Control horizon
mpc_params = [Hp, Hu];           % Hp, Hu
cost_params = [1e-6, 1, 1e-6];    % lambda1, lambda2, epsilon
min_dist = 5;                  % Minimum distance between robots

% Map bounderies
xmin = 0; xmax = 40;
ymin = 0; ymax = 40;
zmin = 0; zmax = 5;
map_bounds = [xmin, xmax, ymin, ymax];

% True initial process states
M_ref = 200;
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
M_0 = 200;
beta_0 = 0.05;
xs_0 = 20;
ys_0 = 20;
M_dot_0 = 0;
beta_dot_0 = 0;
xs_dot_0 = 0;
ys_dot_0 = 0;
z_est_0 = [M_0; M_dot_0; beta_0; beta_dot_0; xs_0; xs_dot_0; ys_0; ys_dot_0]; % Initial guessed state vector

% Initial error covariance matrix
P_0 = zeros(Nx_p);
% P_0(1,1) = (M_ref - M_0)^2;
% P_0(3,3) = (beta_ref - beta_0)^2;
% P_0(5,5) = (xs_ref - xs_0)^2;
% P_0(7,7) = (ys_ref - ys_0)^2;
% P_0(1,1) = 1;
% P_0(5,5) = 1;
% P_0(7,7) = 1;
kf_init = [z_est_0, P_0];

% Initial robot positions, Z0 and control inputs, Uprev
x_1 = zeros(2*M,1); % (assuming M=4, so 2*M values)
for m = 1:M
    x_1(2*m-1) = 1;
    x_1(2*m) = 1;
end
U_prev = zeros(2*M,1);           % (2*M values again)

%% Dynamic Parameter Model
tau_beta = 0.95;       %spread increase parameter (exponential decay of beta)

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
sigma_M = 2; sigma_beta = 0.00001; sigma_x = 0.1; sigma_y = 0.1;
v = [sigma_M sigma_beta sigma_x sigma_y]';
mu_w = zeros(length(v),1);
Q = (inv(T)*G)*diag([sigma_M^2 sigma_beta^2 sigma_x^2 sigma_y^2])*(inv(T)*G)'; % This is the process noise covariance matrix

% Measurement noise
sigma_measurement = 0.1;
mu_r = zeros(M,1);
R = sigma_measurement^2*eye(M); % measurement noise covariance

% % Measurement model - h_i(z,x,y) = (M*beta/pi)*exp(-beta*((x-xs).^2 + (y-ys).^2))
% h_p = @(states,px,py) (states(1)*states(3)/pi)*exp(-states(3)*((px-states(5)).^2 + (py-states(7)).^2));
% h_vec = @(z,x,M) arrayfun(@(m) h_p(z, x(2*m-1), x(2*m)), (1:M)');
% % TESTING
% % h(z_0,x_1)
% % h = zeros(M,1);
% % for m = 1:M
% %     h(m) = h_p(z_0, x_1(2*m-1), x_1(2*m));
% % end
% % h
% 
% % Functions used in computing the Jacobian
% H_m1 = @(z, px, py) h_p(z, px, py)/z(1);
% H_m3 = @(z, px, py) h_p(z, px, py)*(1/z(3) - ((px - z(5))^2 + (py - z(7))^2));
% H_m5 = @(z, px, py) h_p(z, px, py)*(-2*z(3)*(z(5) - px));
% H_m7 = @(z, px, py) h_p(z, px, py)*(-2*z(3)*(z(7) - py));
% H_m = @(z, px, py) [H_m1(z, px, py), 0, H_m3(z, px, py), 0, H_m5(z, px, py), 0, H_m7(z, px, py), 0];
% 
% % Dynamically construct a function that computes the Jacobian for a
% % specific state and robot position
% H = @(z, x, M) cell2mat(...
%             arrayfun(...
%             @(m) H_m(z, x(2*m-1), x(2*m)), (1:M)', 'UniformOutput', false));

% %TESTING
% H(z_0,x_1,M)
%get_H_jaco(z,x_1,M,sc)
% H = zeros(M,Nx_p);
% for m=1:M
%     x_pos = x_1(2*m-1);
%     y_pos = x_1(2*m);
%     y = h_p(z_0, x_pos, y_pos);
%     H(m,1) = y/z_0(1);
%     H(m,3) = -y*((x_pos-z_0(5))^2 + (y_pos-z_0(7))^2) + y/z_0(3);
%     H(m,5) = -2*y*z_0(3)*(z_0(5)-x_pos);
%     H(m,7) = -2*y*z_0(3)*(z_0(7)-y_pos);
% end
% H

% Transformation for scaling
% syms k_M k_Mdot k_beta k_betadot k_xs k_xsdot k_ys k_ysdot real
% T = [k_M 0 0 0 0 0 0 0;
%      0 k_Mdot 0 0 0 0 0 0;
%      0 0 k_beta 0 0 0 0 0;
%      0 0 0 k_betadot 0 0 0 0;
%      0 0 0 0 k_xs 0 0 0;
%      0 0 0 0 0 k_xsdot 0 0;
%      0 0 0 0 0 0 k_ys 0;
%      0 0 0 0 0 0 0 k_ysdot];

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
e(:,1) = ones(M,1)*0.5;     % Robots fully charged at time 0
charge_control = ones(M,K+1);
charge_control(:,1) = zeros(M,1);

% Energy dynamics
t1 = 0.5;
o1 = 5;
charge_rate = 0.2;
en_charge = @(x,y) (1/(1+exp(t1*(x-o1))))*(1/(1+exp(t1*(y-o1))))*charge_rate; % Energy charging function
en_cons = @(v) 0.0001*v^2 + 0.05; % Energy consumption function

%power = @(x,y,v,ctrl) ctrl*(1/(1+exp(t1*(x-o1))))*(1/(1+exp(t1*(y-o1))))*charge_rate - (0.0001*v^2 + 0.05);
%power = @(v) 0.0001*v^2 + 0.05;
%power = 0.05;


% Parameters
charger_x = 0;
charger_y = 0;
charger_r = 5;

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

get_H_jaco(z(:,1),[21;21;21;21;21;21],M,sc)

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
for m = 1:M
    robot_true(m) = plot(ax1, x_1(2*m-1), x_1(2*m), 'o', ...
                                    'MarkerSize', 10, ...
                                    'MarkerFaceColor', colors{m}, ...
                                    'MarkerEdgeColor', 'k', ...
                                    'LineWidth', 1);
    robot_est(m) = plot(ax2, x_1(2*m-1), x_1(2*m), 'o', ...
                                    'MarkerSize', 10, ...
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

z_real = zeros(Nx_p,K+1);
z_est_real = zeros(Nx_p,K+1);

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
    EN_params = {e(:,k), en_charge, en_cons, charger_x, charger_y, charger_r};
    
    % Compute optimal stuff
    [X_opt, U_opt, P_trace, sol_prev, charge_control(:,k)] = MPC_func(ROB_params, KF_params, EN_params, mpc_params, cost_params, M, map_bounds, min_dist, sim_params, 2, do_warm_start, sol_prev, sc);
    u_opt = U_opt(:,2);
    do_warm_start = 1; % set warm start to 1 after first iteration

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