%Radiation/chemical leak (running example) model parameter estimation
clc; clear; close all;
%% Adjustable Parameters
% General simulation parameters
M = 3;                  % Number of robots
dt = 0.1;               % Sampling period [s]
sim_time = 100;          % Simulation time [s]
K = sim_time/dt;        % Total # of simulation steps
Ts = 0.8;                 % MPC sampling period
sim_params = [Ts, dt];
state_plot_selec = [5,7];            % Select states to plot
error_cv_selec = [1,3,5,7];

% MPC parameters
Hp = 6;             % Prediction horizon
Hu = 3;             % Control horizon
mpc_params = [Hp, Hu];           % Hp, Hu
cost_params = [1e-6, 1, 1e-6];    % lambda1, lambda2, epsilon
min_dist = 0.5;                  % Minimum distance between robots

% Map bounderies
xmin = 0; xmax = 40;
ymin = 0; ymax = 40;
zmin = 0; zmax = 20;
map_bounds = [xmin, xmax, ymin, ymax];

% True initial process states
M_ref = 250;
M_dot_ref = 0;
beta_ref = 0.05;
beta_dot_ref = 0;
xs_ref = 15;
xs_dot_ref = 0;
ys_ref = 25;
ys_dot_ref = 0;
z_0 = [M_ref; M_dot_ref; beta_ref; beta_dot_ref; xs_ref; xs_dot_ref; ys_ref; ys_dot_ref]; % Initial true state vector
Nx_p = length(z_0); % Number of states in process state vector

% Guessed initial process states
M_0 = 100;
M_dot_0 = 0;
beta_0 = 0.05;
beta_dot_0 = 0;
xs_0 = 10;
xs_dot_0 = 0;
ys_0 = 10;
ys_dot_0 = 0;
z_est_0 = [M_0; M_dot_0; beta_0; beta_dot_0; xs_0; xs_dot_0; ys_0; ys_dot_0]; % Initial guessed state vector

% Initial error covariance matrix
P_0 = zeros(Nx_p);
P_0(1,1) = (M_ref - M_0)^2;
P_0(3,3) = (beta_ref - beta_0)^2;
P_0(5,5) = (xs_ref - xs_0)^2;
P_0(7,7) = (ys_ref - ys_0)^2;
% P_0(1,1) = 1;
% P_0(5,5) = 1;
% P_0(7,7) = 1;
kf_init = [z_est_0, P_0];

% Initial robot positions, Z0 and control inputs, Uprev
x_1 = zeros(2*M,1); % (assuming M=4, so 2*M values)
for m = 1:M
    x_1(2*m-1) = xs_0;
    x_1(2*m) = ys_0;
end
U_prev = zeros(2*M,1);           % (2*M values again)

%% Dynamic Parameter Model
tau_beta = 0.999;       %spread increase parameter (exponential decay of beta)

% Inputs
Nu_p = 2;
beta0 = 0.05; %steady-state value for beta
u_beta = beta0*(1-tau_beta)/dt; % constant input to achieve beta0 in steady-state
u_M = 10; % Leakage rate offset
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
sigma_M = 1; sigma_beta = 0.00001; sigma_x = 0.1; sigma_y = 0.1;
v = [sigma_M sigma_beta sigma_x sigma_y]';
mu_w = zeros(length(v),1);
Q = G*diag([sigma_M^2 sigma_beta^2 sigma_x^2 sigma_y^2])*G'; % This is the process noise covariance matrix

% Measurement noise
sigma_measurement = 0.001;
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
% get_H_jaco(z_0,x_1,M)
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

% Everyting recomputed for the MPC
A_MPC = A_func(Ts);
B_MPC = B_func(Ts);
Q_MPC = G_func(Ts)*diag([sigma_M^2 sigma_beta^2 sigma_x^2 sigma_y^2])*G_func(Ts)';

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

%% Create Required Vectors and Matrices
% MPC
x = zeros(2*M,K+1);         % Robot positions

% KALMAN FILTER           
z = zeros(Nx_p,K+1);           % True process states (Nx_p)
z_est = zeros(Nx_p,K+1);       % Estimated process states by Kalman filter
P = zeros(Nx_p,Nx_p,K+1);         % Error covariance matrix
y = zeros(M,K+1);           % Measurements of process

% Put initial conditions into vectors
z(:,1) = z_0;
z_est(:,1) = z_est_0;
x(:,1) = x_1;
x(:,2) = x_1;
P(:,:,1) = P_0;

%% Plotting preparation
% Initializations
x_axis = linspace(xmin,xmax,50);    % Process x-axis
y_axis = linspace(ymin,ymax,50)';   % Process y-axis
z_values = get_h(z_0,x_axis,y_axis);    % Process values in defined area
t_vec = (0:K) * dt;               % Time vector for plotting
colors = {'b','r','g','k','y','c','m','b'};     % valid MATLAB colors

% Initialize plotting objects that can be updated 
process_plot = gobjects(1,1);         

% Initial process plot setup
figure(1)
set(gcf, 'Position',  [100, 200, 1200, 500])
subplot(1,2,1)
hold on;
axis([xmin xmax ymin ymax zmin zmax]);
xlabel('X'); ylabel('Y');
process_plot(1) = mesh(x_axis,y_axis,z_values, 'FaceAlpha','0.7', 'EdgeAlpha','0.7');
view([45 45])

% Plot robots on process plot (with their inital conditions)
robot_quivers = gobjects(M,1);
for m = 1:M
    robot_quivers(m) = quiver(x_1(2*m-1), x_1(2*m),...
        0.6, 0.6,...
        'Color', colors{m}, 'LineWidth', 4, 'MaxHeadSize', 4);
end

% Plot selected states
z_selected = zeros(length(state_plot_selec),K+1);          % Vector for storage of selected true process states
z_est_selected = zeros(length(state_plot_selec),K+1);      % Vector for storage of selected estimated process states
z_selected_plot = gobjects(size(z_selected,1),1);
z_est_selected_plot = gobjects(size(z_selected,1),1);
subplot(1,2,2)
hold on
for i = 1:length(state_plot_selec)
    s = state_plot_selec(i);
    z_selected(i,1) = z_0(s);   % Plug in initial states as first column
    z_est_selected(i,1) = z_est_0(s);   % Plug in initial guessed states as first column
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

for k=2:K+1
    clc
    % Update true process
    w = G*normrnd(mu_w,v);
    z(:,k) = A*z(:,k-1) + B*u(:,k-1) + w;
    
    % Update robot positions
    for m = 1:M
        % Update robot positions
        % x(m,k) = z(m,k-1) + normrnd(0, 1);
        % x(2*m-1,k) = z(5,k-1) + normrnd(0, 2);
        % x(2*m,k) = z(7,k-1) + normrnd(0, 2);
        %Randomly update speed/turn rate
        if rand() < 0.5
            robot_v(m) = max_speed * rand();
            robot_omega(m) = (rand()-0.5)*2*max_omega;
        end

        % Update orientation
        robot_theta(m) = mod(robot_theta(m) + robot_omega(m)*dt, 2*pi);

        % Proposed new position
        x(2*m-1,k) = x(2*m-1,k-1) + robot_v(m)*cos(robot_theta(m))*dt;
        x(2*m,k) = x(2*m,k-1) + robot_v(m)*sin(robot_theta(m))*dt;
        % TODO (right now we cheated using MPC output)
    end
    %x(:,k) = A_r*x(:,k) + B_r*u_opt;

    % Take measurements
    r = mvnrnd(mu_r',R)';
    y(:,k) = get_h_vec(z(:,k), x(:,k), M) + r;
    
    % Run Kalman filter iteration
    [z_est(:,k), P(:,:,k), z_hat(:,k)] = EKF(z_est(:,k-1), P(:,:,k-1), A, B, u(:,k-1), Q, y(:,k), R, x(:,k));
   
    % Package all Kalman filter information
    u_mpc = u(:,k:k+Hp);                % This could be avoided by removing DeltaT in B and in the u_beta calc.
    u_mpc(2,:) = u_mpc(2,:)*(dt/Ts);    % Scale input to match MPC sampling period
    KF_params = {z_est(:,k), P(:,:,k), A_MPC, B_MPC, u_mpc, Q_MPC, R, error_cv_selec};

    % Package all robot information
    ROB_params = {x(:,k), A_r_MPC, B_r_MPC, u_opt};
    
    % Compute optimal stuff
    %kf_init = [z_est(:,k), P(:,:,k)];
    %rob_init = [x(:,k), u_opt];
    [X_opt, U_opt, P_trace] = MPC_func(ROB_params, KF_params, mpc_params, cost_params, M, map_bounds, min_dist, sim_params, 3);
    u_opt = U_opt(:,2);
    x(:,k+1) = X_opt(:,2);

    % Update process plot
    subplot(1,2,1)
    z_values = get_h(z(:,k),x_axis,y_axis);
    set(process_plot(1), 'XData', x_axis, 'YData', y_axis, 'ZData', z_values)
    
    % Update robot positions
    for m = 1:M
        set(robot_quivers(m),'XData', x(2*m-1,k),'YData', x(2*m,k));
    end
    
    %update state and estimate plot
    for i = 1:length(state_plot_selec)
        s = state_plot_selec(i);
        z_selected(i,k) = z(s,k);
        z_est_selected(i,k) = z_est(s,k);
        set(z_selected_plot(i), 'XData', t_vec(1:k),'YData', z_selected(i,1:k))
        set(z_est_selected_plot(i), 'XData', t_vec(1:k), 'YData', z_est_selected(i,1:k))
    end
    drawnow limitrate
end

%% Plots 

% OLD STUFF
    % set(z_est_selected_plot(2), 'XData', t_vec(1:k), 'YData', z_est(3,1:k)) % Beta
    % set(z_est_selected_plot(3), 'XData', t_vec(1:k), 'YData', z_est(5,1:k)) % X position
    % set(z_est_selected_plot(4), 'XData', t_vec(1:k), 'YData', z_est(7,1:k)) % Y position

    % % Measurement update
    % for i = 1:M
    %     % Update robot positions
    %     x_pos = x_opt(2*i-1);
    %     y_pos = x_opt(2*i);
    %     pos_store(k,:,i) = [x_pos, y_pos];
    % 
    %     % Take measurements
    %     v = normrnd(0,sigma_measurement);
    %     measurements(i,k) = h(z(:,k), x_pos, y_pos) + v;
    % 
    %     % Make measurement predictions
    %     y_pred(i) = h(z_pred, x_pos, y_pos);
    % 
    %     % Linearize
    %     H(i,1) = y_pred(i)/z_pred(1);
    %     H(i,3) = -y_pred(i)*((x_pos-z_pred(5))^2 + (y_pos-z_pred(7))^2) + y_pred(i)/z_pred(3);
    %     H(i,5) = -2*y_pred(i)*z_pred(3)*(z_pred(5)-x_pos);
    %     H(i,7) = -2*y_pred(i)*z_pred(3)*(z_pred(7)-y_pos);
    % end
    % y_error = measurements(:,k) - y_pred;
    % K = P_pred*H'*(H*P_pred*H' + R)^-1;
    % 
    % estimates(:,k) = z_pred + K*y_error;
    % P_update = (eye(size(z,1)) - K*H)*P_pred*(eye(size(z,1)) - K*H)' + K*R*K';
    % 
    % % Time update
    % z_pred = A*estimates(:,k) + B*u(:,k);
    % P_pred = A*P_update*A' + Q;
    % 
    % kf_init = [z_pred, P_pred];
    % rob_init = [x_opt, u_opt];
    % [X_opt, U_opt, P_trace] = MPC_func(rob_init, kf_init, mpc_params, cost_params, M, map_bounds, min_dist, sim_params, 3);
    % u_opt = U_opt(:,2);
    % x_opt = X_opt(:,2);

% Intensity distribution animation
% nx_vis = 50; ny_vis = 50;
% x_vis = linspace(xmin,xmax,nx_vis);
% y_vis = linspace(ymin,ymax,ny_vis);
% [X_vis, Y_vis] = meshgrid(x_vis, y_vis);

% figure;
% for i = 1:K
%     I_max_true = max(z(1,i))*max(z(3,i))/pi;
%     I_max_est = max(estimates(1,i))*max(estimates(3,i))/pi;
%     I_max = max(I_max_true, I_max_est);
% 
%     I_min_true = min(z(1,i))*min(z(3,i))/pi;
%     I_min_est = min(estimates(1,i))*min(estimates(3,i))/pi;
%     I_min = min(I_min_true, I_min_est);
%     if I_min >= 0
%         I_min = 0;
%     end
% 
%     q_true = (z(1,i)*z(3,i)/pi)*exp(-z(3,i)*((X_vis - z(5,i)).^2 + (Y_vis - z(7,i)).^2));
%     q_est = (estimates(1,i)*estimates(3,i)/pi)*exp(-estimates(3,i)*((X_vis - estimates(5,i)).^2 + (Y_vis - estimates(7,i)).^2));
% 
%     subplot(1,2,1);
%     cla;
%     hold on;
%     imagesc(x_vis, y_vis, q_true);
%     for j = 1:M
%         pos = pos_store(i,:,j); % drone position at step i
%         plot(pos(1), pos(2), 'ko', 'MarkerFaceColor','w', 'MarkerSize',10);
%     end
%     hold off;
%     set(gca,'YDir','normal');
%     %caxis([I_min I_max]); 
%     colorbar;
%     % title('True Radiation Field');
%     title(sprintf('True Radiation Field @t=%.3f',i*dt));
%     xlabel('x'); xlim([xmin, xmax]);
%     ylabel('y'); ylim([ymin, ymax]);
% 
%     subplot(1,2,2);
%     cla;
%     hold on;
%     imagesc(x_vis, y_vis, q_est);
%     for j = 1:M
%         pos = pos_store(i,:,j); % drone position at step i
%         plot(pos(1), pos(2), 'ko', 'MarkerFaceColor','w', 'MarkerSize',10);
%     end
%     hold off;
%     set(gca,'YDir','normal');
%     %caxis([I_min I_max]); 
%     colorbar;
%     % title('Estimated Radiation Field');
%     title(sprintf('Estimated Radiation Field @t=%.3f',i*dt));
%     xlabel('x'); xlim([xmin, xmax]);
%     ylabel('y'); ylim([ymin, ymax]);
% 
%     drawnow limitrate;
%     pause((1/3)*dt);
% end

% Plot state estimates
% figure;
% for i = 1:size(z,1)
%     subplot(size(z,1)/2,2,i);
%     hold on;
%     plot(dt*(0:K-1), z(i,1:K), DisplayName=sprintf("z_{%d}", i));
%     plot(dt*(0:K-1), estimates(i,:), '--', DisplayName=sprintf("z_{%d,est}", i));
%     hold off;
%     legend;
%     xlabel("Time [s]");
%     grid on;
% end
% 
% figure;
% hold on;
% for j = 1:M
%     plot(pos_store(:,1,j), pos_store(:,2,j), '-o', DisplayName=sprintf("Robot %d",j));
%     plot(pos_store(1,1,j), pos_store(1,2,j), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 3, 'HandleVisibility', 'off');
%     plot(pos_store(end,1,j), pos_store(end,2,j), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3, 'HandleVisibility', 'off');
% end
% hold off;
% legend;
% xlabel("x");
% ylabel("y");
% xlim([xmin, xmax]);
% ylim([ymin ymax]);
% grid on;


