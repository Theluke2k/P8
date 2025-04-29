%Radiation/chemical leak (running example) model parameter estimation
clc; clear; close all;

%% Process simulation parameters
N = 3;                  % Number of robots
dt = 0.1;               % Time step [s]
sim_time = 20;          % Simulation time [s]
T = sim_time/dt;        % Total # of simulation steps

% Map bounderies
xmax = 20;
ymax = 20;
xmin = 0;
ymin = 0;

% True initial model parameters
M_ref = 200;
beta_ref = 0.1;
xs_ref = 10;
ys_ref = 10;

% Initial derivative states
M_dot = 0;
beta_dot = 0;
xs_dot = 0;
ys_dot = 0;

% Initial true state vector
z = zeros(8,T+1);
z(:,1) = [M_ref; M_dot; beta_ref; beta_dot; xs_ref; xs_dot; ys_ref; ys_dot]; %initial state vector

% Initial model parameters guess
M = 1e-3;
beta = 1e-3;
xs = 8;
ys = 12;

% Initial states and error covariance matrix
z_pred = [M; M_dot; beta; beta_dot; xs; xs_dot; ys; ys_dot]; % Initial prediction for kalman filter
P_pred = 0*eye(size(z,1)); % Initial error covariance

for i = 1:2:3
    P_pred(i,i) = 1e-1;
end

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
sigma_M = 0.2; sigma_beta = 0.0001; sigma_x = 0.01; sigma_y = 0.01;
sigma = [sigma_M sigma_beta sigma_x sigma_y]';
mu_w = zeros(length(sigma),1);
Q = G*diag([sigma_M sigma_beta sigma_x sigma_y])*G';

% Measurement model - h(z,x,y) = (M*beta/pi)*exp(-beta*((x-xs).^2 + (y-ys).^2))
h = @(states,x,y) (states(1)*states(3)/pi)*exp(-states(3)*((x-states(5)).^2 + (y-states(7)).^2));

% Measurement noise
sigma_measurement = 0.1;
R = sigma_measurement*eye(N); % measurement noise covariance

%% Predifined robot trajectories
pos_store = zeros(T, 2, N);

% Robot 1 - btm left to top right
pos_store(:,:,1) = [linspace(xmin,xmax,T)', linspace(ymin,ymax,T)'];

% Robot 2 - top left to btm right
pos_store(:,:,2) = [linspace(xmin,xmax,T)', linspace(ymax,ymin,T)'];

% Robot 3 - Circle around the center
theta = linspace(0,2*pi,T)';
r = 1.5;
center = [10,10];
pos_store(:,:,3) = [center(1) + r*cos(theta), center(2) + r*sin(theta)];

%% Storage and plotting for measurements and kalman filter
measurements = zeros(N,T);
estimates = zeros(size(z,1),T);
y_pred = zeros(N,1); % Storage for output prediction

%% Simulation loop
H = zeros(N,size(z,1)); % Initialize linearized H matrix
for n=1:T
    % Simulate true process
    w = G*normrnd(mu_w,sigma);
    z(:,n+1) = A*z(:,n) + B*u(:,n) + w;

    % Measurement update
    for i = 1:N
        % Update robot positions
        x_pos = pos_store(n,1,i);
        y_pos = pos_store(n,2,i);

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
    caxis([I_min I_max]); 
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
    caxis([I_min I_max]); 
    colorbar;
    % title('Estimated Radiation Field');
    title(sprintf('Estimated Radiation Field @t=%.3f',i*dt));
    xlabel('x'); xlim([xmin, xmax]);
    ylabel('y'); ylim([ymin, ymax]);

    drawnow limitrate;
    pause(0.001);
end

% Plot state estimates
figure;
for i = 1:size(z,1)
    subplot(size(z,1)/2,2,i);
    hold on;
    plot(0:T-1, z(i,1:T), DisplayName=sprintf("z_{%d}", i));
    plot(0:T-1, estimates(i,:), '--', DisplayName=sprintf("z_{%d,est}", i));
    hold off;
    legend;
    xlabel("Time [s]");
    grid on;
end
