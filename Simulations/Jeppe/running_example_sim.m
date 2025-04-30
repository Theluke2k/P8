%Radiation/chemical leak (running example) model parameter estimation
clear all
close all
%% Process simulation parameters
N = 3;                  % Number of robots
dt = 0.1;               % Time step [s]
sim_time = 100;          % Simulation time [s]
T = sim_time/dt;        % Total # of simulation steps

% Map bounderies
xmax = 20;
ymax = 20;
xmin = 0;
ymin = 0;

zmin = 0;
zmax = 20;

%Initial model parameters
M = 250;        % Total leak mass
beta = 0.05;    % Spread
xs = 10;        % Source x-coordinate
ys = 10;        % Source y-coordinate

%Initial derivative states
M_dot = 0;
beta_dot = 0;
xs_dot = 0;
ys_dot = 0;

%% Dynamic parameter model
z = [M; M_dot; beta; beta_dot; xs; xs_dot; ys; ys_dot]; %initial state vector
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

B = zeros(length(z),2);
B(1,1) = dt; B(3,2) = dt;

%Process noise
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

%% Measurement model
sigma_measurement = 0.1;
R = sigma_measurement*eye(N); %measurement noise covariance
%q_measurement = @(x,y) (M*beta/pi)*exp(-beta*((x-xs).^2 + (y-ys).^2));

%% Robot Initialization
% Generate initial positions in a diagonal line starting from (10, 10)
%robot_pos = [10, 10] + (0:M-1)'*[1, 1]; % MOHAMAD
robot_posX = (xmax - xmin)*rand(N,1) + xmin;
robot_posY = (ymax - ymin)*rand(N,1) + ymin;
robot_pos = [robot_posX robot_posY];

% Set initial orientations to be equally spaced around a circle
%robot_theta = (0:M-1)'*(2*pi/M); % MOHAMAD
robot_theta = rand(N,1) * 2*pi;
robot_v = zeros(N,1);          % Linear velocities
robot_omega = zeros(N,1);      % Angular velocities
%% Movement Parameters
max_speed = 0.8*2;        % Maximum robot speed
max_omega = pi/2;       % Maximum angular velocity
change_prob = 0.1;      % Probability of changing direction

%% Initialize process plotting
q = @(x,y) (M*beta/pi)*exp(-beta*((x-xs).^2 + (y-ys).^2));
x_plot = linspace(xmin,xmax,50);
y_plot = linspace(ymin,ymax,50)';
q_plot = q(x_plot,y_plot);

t_vec = (0:T-1) * dt;         % Time vector for plotting
parameters_time = zeros(T,4); % Vector for storage of parameters (M, beta, xs, ys)
amplitude_time = zeros(T,1);

field_plot = gobjects(1,1);
parameter_plot = gobjects(4,1);
amplitude_plot = gobjects(1,1);

%% Storage and plotting for measurements and kalman filter
measurements = zeros(N,T);
estimates = zeros(length(z),T);

z_pred = z; % Initial prediction for kalman filter
P_pred = zeros(length(z)); % Initial error covariance
y_pred = zeros(N,1); %storage for output prediction
H = zeros(N,length(z));

% Struct for UKF parameters (default values from TK's slides)
UKFparams.alpha = 1;
UKFparams.kappa = 2;
UKFparams.beta = 0;
UKFparams.n = length(z);
UKFparams.lambda = UKFparams.alpha^2*(UKFparams.n + UKFparams.kappa) + UKFparams.n;
UKFparams.k = sqrt(UKFparams.n + UKFparams.lambda);

estimate_plot = gobjects(4,1);

%% Generate plots

figure(1)
set(gcf, 'Position',  [100, 200, 1200, 500])

subplot(1,2,1)
hold on;
axis([xmin xmax ymin ymax zmin zmax]);
xlabel('X'); ylabel('Y');
field_plot(1) = mesh(x_plot,y_plot,q_plot, 'FaceAlpha','0.7', 'EdgeAlpha','0.7');
view([45 45])

robot_quivers = gobjects(N,1);
for i = 1:N
    robot_quivers(i) = quiver(robot_pos(i,1), robot_pos(i,2),...
        0.5*cos(robot_theta(i)), 0.5*sin(robot_theta(i)),...
        'r', 'LineWidth', 4, 'MaxHeadSize', 4);
end

subplot(1,2,2)
parameter_plot(1) = plot(nan,nan,'DisplayName','$\dot M$','Color','b');
hold on
parameter_plot(2) = plot(nan,nan,'DisplayName','$\beta$','Color','r');
parameter_plot(3) = plot(nan,nan,'DisplayName','$x_s$','Color','g');
parameter_plot(4) = plot(nan,nan,'DisplayName','$y_s$','Color','k');
%amplitude_plot(1) = plot(nan,nan);
estimate_plot(1) = plot(nan,nan,'DisplayName','$\hat{\dot M}$','Color','b','LineStyle','--');
estimate_plot(2) = plot(nan,nan,'DisplayName','$\hat{\beta}$','Color','r','LineStyle','--');
estimate_plot(3) = plot(nan,nan,'DisplayName','$\hat{x_s}$','Color','g','LineStyle','--');
estimate_plot(4) = plot(nan,nan,'DisplayName','$\hat{y_s}$','Color','k','LineStyle','--');
xlabel('Time [s]')
legend('Interpreter','latex','Location','west')

%% Simulation loop
for n=1:T

    %simulate process
    w = G*normrnd(mu_w,sigma);
    z = A*z + B*u(:,n) + w;
    parameters_time(n,1) = z(2) + u(1,n); % Store leak rate (M_dot)
    parameters_time(n,2) = z(3); % Store beta
    parameters_time(n,3) = z(5); % Store xs
    parameters_time(n,4) = z(7); % Store ys
    amplitude_time(n,1) = z(1)*z(3)/pi;
    q = @(x,y) (z(1)*z(3)/pi)*exp(-z(3)*((x-z(5)).^2 + (y-z(7)).^2));

    % 2) Update robot positions (with boundary check)
    for i = 1:N
        % Randomly update speed/turn rate
        if rand() < change_prob
            robot_v(i) = max_speed * rand();
            robot_omega(i) = (rand()-0.5)*2*max_omega;
        end
        
        % Update orientation
        robot_theta(i) = mod(robot_theta(i) + robot_omega(i)*dt, 2*pi);
        
        % Proposed new position
        new_x = robot_pos(i,1) + robot_v(i)*cos(robot_theta(i))*dt;
        new_y = robot_pos(i,2) + robot_v(i)*sin(robot_theta(i))*dt;
        
        % If out of bounds, freeze movement
        if (new_x < xmin || new_x > xmax || new_y < ymin || new_y > ymax)
            robot_v(i) = 0;  % Robot "freezes"
        else
            % Update position only if within boundary
            robot_pos(i,1) = new_x;
            robot_pos(i,2) = new_y;
        end
        %Take measurement
        v = normrnd(0,sigma_measurement);
        measurements(i,n) = q(robot_pos(i,1), robot_pos(i,2)) + v;
    end

    % Excute EKF

    q_pred = @(x,y) (z_pred(1)*z_pred(3)/pi)*exp(-z_pred(3)*((x-z_pred(5)).^2 + (y-z_pred(7)).^2));
    
    for i=1:N
        x_pos = robot_pos(i,1);                 % For readability
        y_pos = robot_pos(i,2);                 % For readability
        y_pred(i) = q_pred(x_pos, y_pos);       % Predicted measurement
        H(i,1) = y_pred(i)/z_pred(1);           % dq/dM (i,1) in Jacobian
        H(i,3) = -y_pred(i)*((x_pos-z_pred(5))^2 + (y_pos-z_pred(7))^2) + y_pred(i)/z_pred(3);      % dq/dM (i,1) in Jacobian
        H(i,5) = -2*y_pred(i)*z_pred(3)*(z_pred(5)-x_pos);
        H(i,7) = -2*y_pred(i)*z_pred(3)*(z_pred(7)-y_pos);
    end
    y_error = measurements(:,n) - y_pred;
    K = P_pred*H'*(H*P_pred*H' + R)^-1;
    estimates(:,n) = z_pred + K*y_error;
    P_update = (eye(length(z)) - K*H)*P_pred*(eye(length(z)) - K*H)' + K*R*K';

    %Time update step
    z_pred = A*estimates(:,n) + B*u(:,n);
    P_pred = A*P_update*A' + Q;

    %Update plots
    subplot(1,2,1)
    %Update field plot

    q_plot = q(x_plot,y_plot);
    set(field_plot(1), 'XData', x_plot, 'YData', y_plot, 'ZData', q_plot)
    
    %Update robot animation
    for i = 1:N
    set(robot_quivers(i), ...
        'XData', robot_pos(i,1), ...
        'YData', robot_pos(i,2), ...
        'UData', 0.5*cos(robot_theta(i)), ...
        'VData', 0.5*sin(robot_theta(i)));
    end
    
    subplot(1,2,2)
    %update parameter and estimates plot
    for i=1:4
    set(parameter_plot(i), 'XData', t_vec(1:n),'YData', parameters_time(1:n,i))
    end
    %set(amplitude_plot(1), 'XData', t_vec(1:n),'YData', amplitude_time(1:n))
    set(estimate_plot(1), 'XData', t_vec(1:n), 'YData', estimates(2,1:n)+u(1,1:n)) %leak rate M_dot
    set(estimate_plot(2), 'XData', t_vec(1:n), 'YData', estimates(3,1:n)) %beta
    set(estimate_plot(3), 'XData', t_vec(1:n), 'YData', estimates(5,1:n)) %x position
    set(estimate_plot(4), 'XData', t_vec(1:n), 'YData', estimates(7,1:n)) %y position
    
    % drawnow limitrate;
    pause(0.02);
end

%% Make new figures of parameters and estimates
figure(2)
plot(t_vec,parameters_time(:,1))
hold on
plot(t_vec,estimates(2,:)+u(1,:))
xlabel('Time [s]')
legend('Leak rate', 'Leak rate estimate', Location='best')

figure(3)
plot(t_vec,parameters_time(:,2))
hold on
plot(t_vec,estimates(3,:))
xlabel('Time [s]')
legend('Beta', 'Beta estimate', Location='best')

figure(4)
plot(t_vec,parameters_time(:,3))
hold on
plot(t_vec,estimates(5,:))
xlabel('Time [s]')
legend('x_s', 'x_s estimate', Location='best')

figure(5)
plot(t_vec,parameters_time(:,4))
hold on
plot(t_vec,estimates(7,:))
xlabel('Time [s]')
legend('y_s', 'y_s estimate', Location='best')
