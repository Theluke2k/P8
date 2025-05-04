%Radiation/chemical leak (running example) model parameter estimation
clc; clear; close all;
addpath('C:\Program Files\MATLAB/casadi-3.7.0-windows64-matlab2018b')
%% Process simulation parameters
N = 3;                  % Number of robots
dt = 0.5;               % Time step [s]
sim_time = 10;          % Simulation time [s]
T = sim_time/dt;        % Total # of simulation steps
Ts = 1;                 % MPC sampling time
sim_params = [Ts, dt];

% MPC and cost params
Hp = 6; Hu = 3;
mpc_params = [Hp, Hu];           % Hp, Hu
cost_params = [1e-6, 1, 1e-6];    % lambda1, lambda2, epsilon
min_dist = 0.5;                  % Minimum distance between robots

% Map bounderies
xmin = -40; xmax = 40;
ymin = -40; ymax = 40;
map_bounds = [xmin, xmax, ymin, ymax];

% True initial model parameters
M_ref = 200;
beta_ref = 0.1;
xs_ref = 20;
ys_ref = 20;

% Initial derivative states
M_dot = 0;
beta_dot = 0;
xs_dot = 20;
ys_dot = 20;

% Initial true state vector
z = zeros(8,T+1);
z(:,1) = [M_ref; M_dot; beta_ref; beta_dot; xs_ref; xs_dot; ys_ref; ys_dot]; %initial state vector

% Initial model parameters guess
M = 200;
beta = 0.1;
xs = 20;
ys = 20;

% Initial states and error covariance matrix
z_pred = [M; M_dot; beta; beta_dot; xs; xs_dot; ys; ys_dot]; % Initial prediction for kalman filter

% Initial error covariance
P_pred = zeros(size(z,1));
P_pred(1,1) = (M_ref - M)^2;
P_pred(3,3) = (beta_ref - beta)^2;
P_pred(5,5) = (xs_ref - xs)^2;
P_pred(7,7) = (ys_ref - ys)^2;
kf_init = [z_pred, P_pred];

% Initial robot positions, Z0 and control inputs, Uprev
X0 = zeros(2*N,1); % (assuming N=4, so 2*N values)
for j = 1:N
    X0(2*j-1) = 0 + j;
    X0(2*j) = 0;
end
Uprev = zeros(2*N,1);           % (2*N values again)


%% Dynamic parameter model
tau_beta = 0.998;       %spread increase parameter (exponential decay of beta)

% Inputs
beta0 = 0; %steady-state value for beta
u_beta = beta0*(1-tau_beta)/dt; % constant input to achieve beta0 in steady-state
u_M = 2; % Leakage rate offset
u = ones(2,T);
u(1,:) = u(1,:)*u_M;
u(2,:) = u(2,:)*u_beta;

A = [1 dt 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0;
     0 0 tau_beta dt 0 0 0 0;
     0 0 0 1 0 0 0 0;
     0 0 0 0 1 dt 0 0;
     0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 dt;
     0 0 0 0 0 0 0 1];

B = zeros(size(z,1),2);
B(1,1) = dt; B(3,2) = dt;

% Process noise
G = [0.5*dt^2 0 0 0;
    dt 0 0 0;
    0 0.5*dt^2 0 0;
    0 dt 0 0;
    0 0 0.5*dt^2 0;
    0 0 dt 0;
    0 0 0 0.5*dt^2;
    0 0 0 dt];
sigma_M = 0.5; sigma_beta = 0.001; sigma_x = 0.2; sigma_y = 0.2;
sigma = [sigma_M sigma_beta sigma_x sigma_y]';
mu_w = zeros(length(sigma),1);
Q = G*diag([sigma_M sigma_beta sigma_x sigma_y])*G';

% Measurement model - h(z,x,y) = (M*beta/pi)*exp(-beta*((x-xs).^2 + (y-ys).^2))
h = @(states,x,y) (states(1)*states(3)/pi)*exp(-states(3)*((x-states(5)).^2 + (y-states(7)).^2));

% Measurement noise
sigma_measurement = 0.1;
R = sigma_measurement*eye(N); % measurement noise covariance

%% Predifined robot trajectories
% pos_store = zeros(T, 2, N);
% 
% % Robot 1 - btm left to top right
% pos_store(:,:,1) = [linspace(xmin,xmax,T)', linspace(ymin,ymax,T)'];
% 
% % Robot 2 - top left to btm right
% pos_store(:,:,2) = [linspace(xmin,xmax,T)', linspace(ymax,ymin,T)'];
% 
% % Robot 3 - Circle around the center
% theta = linspace(0,2*pi,T)';
% r = 1.5;
% center = [10,10];
% pos_store(:,:,3) = [center(1) + r*cos(theta), center(2) + r*sin(theta)];

%% Storage and plotting for measurements and kalman filter
measurements = zeros(N,T);
estimates = zeros(size(z,1),T);
y_pred = zeros(N,1); % Storage for output prediction

pos_store = zeros(T,2,N);

%% Simulation loop
x_opt = X0;
u_opt = Uprev;
% rob_init = [x_opt, u_opt];

H = zeros(N,size(z,1)); % Initialize linearized H matrix
for n=1:T
    disp("")
    % Simulate true process
    w = G*normrnd(mu_w,sigma);
    z(:,n+1) = A*z(:,n) + B*u(:,n) + w;

    % Measurement update
    for i = 1:N
        % Update robot positions
        x_pos = x_opt(2*i-1);
        y_pos = x_opt(2*i);
        pos_store(n,:,i) = [x_pos, y_pos];

        % Take measurements
        v = normrnd(0,sigma_measurement);
        measurements(i,n) = h(z(:,n), x_pos, y_pos) + v;
        
        % Make measurement predictions
        y_pred(i) = h(z_pred, x_pos, y_pos);
        
        % Linearize
        H(i,1) = y_pred(i)/z_pred(1);
        H(i,3) = -y_pred(i)*((x_pos-z_pred(5))^2 + (y_pos-z_pred(7))^2) + y_pred(i)/z_pred(3);
        H(i,5) = -2*y_pred(i)*z_pred(3)*(z_pred(5)-x_pos);
        H(i,7) = -2*y_pred(i)*z_pred(3)*(z_pred(7)-y_pos);
    end
    y_error = measurements(:,n) - y_pred;
    K = P_pred*H'*(H*P_pred*H' + R)^-1;
    
    estimates(:,n) = z_pred + K*y_error;
    P_update = (eye(size(z,1)) - K*H)*P_pred*(eye(size(z,1)) - K*H)' + K*R*K';

    % Time update
    z_pred = A*estimates(:,n) + B*u(:,n);
    P_pred = A*P_update*A' + Q;
    
    kf_init = [z_pred, P_pred];
    rob_init = [x_opt, u_opt];
    [X_opt, U_opt, P_trace] = MPC_func(rob_init, kf_init, mpc_params, cost_params, N, map_bounds, min_dist, sim_params, 3);
    u_opt = U_opt(:,2);
    x_opt = X_opt(:,2);
end

%% Plots 

% Intensity distribution animation
nx_vis = 50; ny_vis = 50;
x_vis = linspace(xmin,xmax,nx_vis);
y_vis = linspace(ymin,ymax,ny_vis);
[X_vis, Y_vis] = meshgrid(x_vis, y_vis);

figure;
for i = 1:T
    I_max_true = max(z(1,i))*max(z(3,i))/pi;
    I_max_est = max(estimates(1,i))*max(estimates(3,i))/pi;
    I_max = max(I_max_true, I_max_est);

    I_min_true = min(z(1,i))*min(z(3,i))/pi;
    I_min_est = min(estimates(1,i))*min(estimates(3,i))/pi;
    I_min = min(I_min_true, I_min_est);
    if I_min >= 0
        I_min = 0;
    end

    q_true = (z(1,i)*z(3,i)/pi)*exp(-z(3,i)*((X_vis - z(5,i)).^2 + (Y_vis - z(7,i)).^2));
    q_est = (estimates(1,i)*estimates(3,i)/pi)*exp(-estimates(3,i)*((X_vis - estimates(5,i)).^2 + (Y_vis - estimates(7,i)).^2));

    subplot(1,2,1);
    cla;
    hold on;
    imagesc(x_vis, y_vis, q_true);
    for j = 1:N
        pos = pos_store(i,:,j); % drone position at step i
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor','w', 'MarkerSize',10);
    end
    hold off;
    set(gca,'YDir','normal');
    %caxis([I_min I_max]); 
    colorbar;
    % title('True Radiation Field');
    title(sprintf('True Radiation Field @t=%.3f',i*dt));
    xlabel('x'); xlim([xmin, xmax]);
    ylabel('y'); ylim([ymin, ymax]);

    subplot(1,2,2);
    cla;
    hold on;
    imagesc(x_vis, y_vis, q_est);
    for j = 1:N
        pos = pos_store(i,:,j); % drone position at step i
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor','w', 'MarkerSize',10);
    end
    hold off;
    set(gca,'YDir','normal');
    %caxis([I_min I_max]); 
    colorbar;
    % title('Estimated Radiation Field');
    title(sprintf('Estimated Radiation Field @t=%.3f',i*dt));
    xlabel('x'); xlim([xmin, xmax]);
    ylabel('y'); ylim([ymin, ymax]);

    drawnow limitrate;
    pause((1/3)*dt);
end

% Plot state estimates
figure;
for i = 1:size(z,1)
    subplot(size(z,1)/2,2,i);
    hold on;
    plot(dt*(0:T-1), z(i,1:T), DisplayName=sprintf("z_{%d}", i));
    plot(dt*(0:T-1), estimates(i,:), '--', DisplayName=sprintf("z_{%d,est}", i));
    hold off;
    legend;
    xlabel("Time [s]");
    grid on;
end

figure;
hold on;
for j = 1:N
    plot(pos_store(:,1,j), pos_store(:,2,j), '-o', DisplayName=sprintf("Robot %d",j));
    plot(pos_store(1,1,j), pos_store(1,2,j), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 3, 'HandleVisibility', 'off');
    plot(pos_store(end,1,j), pos_store(end,2,j), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3, 'HandleVisibility', 'off');
end
hold off;
legend;
xlabel("x");
ylabel("y");
xlim([xmin, xmax]);
ylim([ymin ymax]);
grid on;
