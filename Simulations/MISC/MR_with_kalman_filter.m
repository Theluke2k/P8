clear all
close all

%% Simulation Parameters
M = 5;                  % Number of robots
N = 1;                  % Number of processes
dt = 0.1;               % Time step [s]
sim_time = 60         % Simulation time [s]
show_contour = true;    % Show measurement quality contours
kalman = 1;             % Toggle for using kalman filter (only works for N=1)
show_all_measurements = 0;

% Map bounderies
xmax_m = 20;
ymax_m = 20;
xmin_m = 0;
ymin_m = 0;
sigmoid = @(x,v) 1 ./ (1+exp(-v*x));
%% Process Parameters (Double Integrated Random Walk)
proc_pos = ([xmax_m ymax_m]-[xmin_m ymin_m]).*rand(N,2) + [xmin_m ymin_m];  % Place process positions randomly within the map boundaries.
A = cell(1,N);          % State transition matrices
G = cell(1,N);          % Process noise matrices
C = cell(1,N);          % Measurement matrices
Q = cell(1,N);          % Process noise covariance

for n = 1:N
    A{n} = [1 dt; 0 1];
    G{n} = [dt^2/2; dt];
    C{n} = [1 0];       % Measure position only
    Q{n} = 0.1;         % Process noise variance
end

%% Measurement Noise Function (R increases with distance)
R0 = 0.1;               % Minimum noise variance
R_max = 1000;            % Maximum noise variance
R_func = @(d) 0.01*exp(0.5*d);

%% Robot Initialization
% Generate initial positions in a diagonal line starting from (10, 10)
%robot_pos = [10, 10] + (0:M-1)'*[1, 1]; % MOHAMAD
robot_posX = (xmax_m - xmin_m)*rand(M,1) + xmin_m;
robot_posY = (ymax_m - ymin_m)*rand(M,1) + ymin_m;
robot_pos = [robot_posX robot_posY];

% Set initial orientations to be equally spaced around a circle
%robot_theta = (0:M-1)'*(2*pi/M); % MOHAMAD
robot_theta = rand(M,1) * 2*pi;
robot_v = zeros(M,1);          % Linear velocities
robot_omega = zeros(M,1);      % Angular velocities
%% Movement Parameters
max_speed = 0.8*2;        % Maximum robot speed
max_omega = pi/2;       % Maximum angular velocity
change_prob = 0.1;      % Probability of changing direction

%% Process States (Position and Velocity)
z = zeros(2, N);        % Each column = [position; velocity] for the process

%% Simulation Setup
T = sim_time/dt;        % Number of simulation steps
t_vec = (0:T-1) * dt;   % Time vector for plotting

%% Create storage for time-series data
% true_vals_time(k,n) = true value of process n at time-step k
true_vals_time = zeros(T, N);

% meas_vals_time(k,m,n) = measurement of process n by robot m at time-step k
meas_vals_time = zeros(T, M, N);

%% Initialize kalman filter
if kalman
kalman_estimates = zeros(T,2);
Q_kalman = [Q{1} 0; 0 Q{1}];
C_kalman = zeros(M,2);
for m=1:M
    C_kalman(m,:) = C{1};
end
z0 = [0;0];
P0 = Q_kalman;
z_pred = z0;
P_pred = P0;
R_kalman = zeros(M,M);
y_multi = zeros(M,1);
end
%% Figure with Two Subplots
figure;

% Left subplot: Animation
subplot(1,2,1);
hold on;
axis([xmin_m xmax_m ymin_m ymax_m]);
grid on;
title('Multi-Robot Process Monitoring');
xlabel('X'); ylabel('Y');

% Plot the (fixed) processes
proc_plot = plot(proc_pos(:,1), proc_pos(:,2), 'rx',...
    'MarkerSize', 10, 'LineWidth', 2);

% Initialize robot quiver objects
robot_quivers = gobjects(M,1);
for m = 1:M
    robot_quivers(m) = quiver(robot_pos(m,1), robot_pos(m,2),...
        0.5*cos(robot_theta(m)), 0.5*sin(robot_theta(m)),...
        'b', 'LineWidth', 2, 'MaxHeadSize', 2);
end

% Optionally, add the contour plots for measurement quality
if show_contour
    [X, Y] = meshgrid(xmin_m:0.1:xmax_m, ymin_m:0.1:ymax_m);
    contour_plots = gobjects(N,1);
    colorsProc = lines(N);
    for n = 1:N
        R_inv = 1 ./ arrayfun(@(xx,yy) R_func(norm([xx,yy]-proc_pos(n,:))), X, Y);
        [~, contour_plots(n)] = contour(X, Y, R_inv, 'LineWidth', 1,...
            'EdgeColor', colorsProc(n,:));
    end
end

% Right subplot: Time-series of true processes & measurements
subplot(1,2,2);
hold on; grid on;
title('True Process Values vs. Measurements Over Time');
xlabel('Time [s]');
ylabel('Value');

% We will have N "true process" lines, and M lines for each process's measurement
colors = lines(N);

% Handles for true process lines: one per process
true_val_plots = gobjects(N,1);
% Handles for measurement lines: M for each process => M-by-N array
meas_val_plots = gobjects(M, N);
if kalman
    kalman_plots = gobjects(N,1);
end

for n = 1:N
    % True process line
    true_val_plots(n) = plot(nan, nan, 'Color', colors(n,:),...
        'LineWidth', 2, 'DisplayName', sprintf('Process %d (True)', n));

    % Measurement lines for each robot
    if show_all_measurements
    for m = 1:M
        % Use a dashed style but the same color as the "true" line
        meas_val_plots(m,n) = plot(nan, nan, '--', 'Color', colors(n,:),...
            'LineWidth', 1, 'DisplayName', sprintf('Process %d, Robot %d', n, m));
    end
    end
    if kalman
        kalman_plots(n) = plot(nan, nan, 'Color', 'r',...
        'LineWidth', 2, 'DisplayName', sprintf('Estimate %d', n));
    end
end
legend('Location','best');

%% Simulation Loop
for k = 1:T
    % 1) Update the process states
    for n = 1:N
        w = sqrt(Q{n}) * randn();
        z(:,n) = A{n} * z(:,n) + G{n} * w;
    end
    
    % 2) Update robot positions (with boundary check)
    for m = 1:M
        % Randomly update speed/turn rate
        if rand() < change_prob
            robot_v(m) = max_speed * rand();
            robot_omega(m) = (rand()-0.5)*2*max_omega;
        end
        
        % Update orientation
        robot_theta(m) = mod(robot_theta(m) + robot_omega(m)*dt, 2*pi);
        
        % Proposed new position
        new_x = robot_pos(m,1) + robot_v(m)*cos(robot_theta(m))*dt;
        new_y = robot_pos(m,2) + robot_v(m)*sin(robot_theta(m))*dt;
        
        % If out of bounds, freeze movement
        if (new_x < xmin_m || new_x > xmax_m || new_y < ymin_m || new_y > ymax_m)
            robot_v(m) = 0;  % Robot "freezes"
        else
            % Update position only if within boundary
            robot_pos(m,1) = new_x;
            robot_pos(m,2) = new_y;
        end
    end
    
    % 3) Collect measurements from each robot for each process
    for n = 1:N
        % True process value (position = first element)
        true_vals_time(k,n) = C{n} * z(:,n);
    end
    
    for m = 1:M
        for n = 1:N
            d = norm(robot_pos(m,:) - proc_pos(n,:));
            R = R_func(d);
            noise = sqrt(R) * randn();
            meas_val = C{n} * z(:,n) + noise;
            meas_vals_time(k,m,n) = meas_val;
            if kalman
                R_kalman(m,m) = R;
                y_multi(m) = meas_val;
            end
        end
    end

    % 4) Execute kalman filter
    if kalman
    %Measurement Update
    y_pred = C_kalman*z_pred;
    y_tilde = y_multi - y_pred;
    K_K = P_pred*C_kalman'*(C_kalman*P_pred*C_kalman' + R_kalman)^-1;
    kalman_estimates(k,:) = z_pred + K_K*y_tilde;
    P_update = (eye(2)-K_K*C_kalman)*P_pred*(eye(2)-K_K*C_kalman)' + K_K*R_kalman*K_K';

    %Time update
    x_pred = A{1}*kalman_estimates(k);
    P_pred = A{1}*P_update*A{1}' + Q_kalman;
    end
    
    %% Update plots
    
    % -- (a) Animation on the left subplot
    subplot(1,2,1);
    for m = 1:M
        set(robot_quivers(m), ...
            'XData', robot_pos(m,1), ...
            'YData', robot_pos(m,2), ...
            'UData', 0.5*cos(robot_theta(m)), ...
            'VData', 0.5*sin(robot_theta(m)));
    end
    
    % -- (b) Time-series on the right subplot
    subplot(1,2,2);
    % Update "true" lines
    for n = 1:N
        set(true_val_plots(n), ...
            'XData', t_vec(1:k), ...
            'YData', true_vals_time(1:k,n));
        
        % Update robot measurement lines
        if show_all_measurements
        for m = 1:M
            set(meas_val_plots(m,n), ...
                'XData', t_vec(1:k), ...
                'YData', meas_vals_time(1:k,m,n));
        end
        end
        %Update Kalman filter line
        set(kalman_plots(n), ...
            'XData', t_vec(1:k), ...
            'YData', kalman_estimates(1:k,n));
    end
    
    drawnow limitrate;
    pause(0.01);
end
