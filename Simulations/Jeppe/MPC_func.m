function [X_opt,U_opt,P_trace, sol] = MPC_func(rob_init, kf_init, en_params, mpc_params, cost_params, n_robots, map_bounds, min_dist, sim_params, verbose_opt, do_warm_start, sol_prev)
%% Parameters
import casadi.*

% Map boundaries
xmin = map_bounds(1); xmax = map_bounds(2);
ymin = map_bounds(3); ymax = map_bounds(4);
dmin = min_dist;
vmax = 2;

% MPC parameters
Hp = mpc_params(1); % Prediction horizon
Hu = mpc_params(2); % Control horizon

% MPC cost function parameters
lambda1 = cost_params(1); % Weight for smooth control
lambda2 = cost_params(2); % Weight for uncertainty
eps = cost_params(3); % Regularization parameter

% Simulation parameters
Ts = sim_params(1); % MPC sampling time
dt = sim_params(2); % Simulation sampling time

% MISC
M = n_robots; % Number of robots

%% KALMAN FILTER
% Define parameters to be used
A_p = kf_init{3};     % Process system matrix with MPC dt
B_p = kf_init{4};     % Process input matrix with MPC dt
u_p = kf_init{5};          % Process input
Q_p = kf_init{6};          % Process noise covariance matrix
R_p = kf_init{7};          % Measurement noise covariance matrix
error_cv_selec = kf_init{8};% Selection of diagonal elements to have in cost function
z_est_0 = kf_init{1};    % KF state from k-1
P_0 = kf_init{2};      % KF error covariance matrix fom k-1

% Useful numbers
Nx_p = size(A_p,1);       % Number of states in process model
Nu_p = size(B_p,2);       % Number of inputs in process model

%% Robot model and Tuning
% Tuning for a single robot
% Q = eye(Nx_p);
% %Q(1,1) = (20/500)^2;
% Q(2,2) = 0;
% Q(4,4) = 0;
% Q(6,6) = 0;
% Q(8,8) = 0;
Q_vec = [1;0;1;0;1;0;1;0];
R_single = [0 0;
     0 0];

% Import robot dynamics from parameters
A_r = rob_init{2};    % Robot system matrix
B_r = rob_init{3};    % Robot input matrix

% Save useful numbers
Nx_r = size(A_r,1)/M;   % Number of states for a sigle robot
Nu_r = size(B_r,2)/M;   % Number of inputs for a single robot

% Reformulate to delta u
A_rd = [A_r,B_r; zeros(M*Nu_r, M*Nx_r), eye(M*Nu_r, M*Nu_r)];
B_rd = [B_r; eye(M*Nu_r)];

% Construct full tuning matrices
%Q = kron(eye(Nx_p), Q_single);
R = kron(eye(M), R_single);

% Robot initial conditions
pos0 = rob_init{1};
u0 = rob_init{4};
x0 = [pos0; u0];    % Initial condition for delta u formulation

%% ENERGY
% Extract energy package info
e0 = en_params{1};
power = en_params{2};
charge_rate = en_params{3};
charger_x = en_params{4};
charger_y = en_params{5};
charger_r = en_params{6};

%% MPC object
% Setup opti
opti = casadi.Opti();

% Decision variables (to be optimized over)
du = opti.variable(M*Nu_r, Hu);   % Delta u (change in robot control input)  
x_rd = opti.variable(M*Nx_r + M*Nu_r,Hp+1);            % Robot positions defined as optimization variables but are bounded by constraints  
b_slack = opti.variable(M,Hp+1);
e_low_slack = opti.variable(M,Hp+1);
e_high_slack = opti.variable(M,Hp+1);

% States and varialbes (NOT to be optimized over)
z_est = MX.zeros(Nx_p,Hp+1);          % Kalman filter states
P = MX.zeros(Nx_p,(Hp+1)*Nx_p);       % Kalman filter error covariance matrix
p = MX.zeros(Nx_p,Hp+1);              % KF error covariance diagonal elements


% Enforce initial conditions
opti.subject_to(x_rd(:,1) == x0);
z_est(:,1) = z_est_0;
P(:,1:Nx_p) = P_0;



%% Create Cost Function
% Control input
cost = 0;
for i = 1:Hu
    cost = cost + du(:,i)'*R*du(:,i);
end

% Kalman filter iterations
z_hat_dum = MX.zeros(Nx_p,Hp+1);
K_dum = MX.zeros(Nx_p,(Hp+1)*M);
H_dum = MX.zeros(M,(Hp+1)*Nx_p);
P_hat_dum = MX.zeros(Nx_p,(Hp+1)*Nx_p);

for i = 2:Hp+1
    % % Robot dynamics
    % if i <= Hu+1
    %     x_rd(:,i) = A_rd*x_rd(:,i-1) + B_rd*du(:,i);
    % else
    %     x_rd(:,i) = A_rd*x_rd(:,i-1);
    % end

    % Define temp. variables
    H = MX.zeros(M,Nx_p);
    h_vec = MX.zeros(M,1);
    K = MX.zeros(Nx_p,M);
    z_hat = MX.zeros(Nx_p,1);
    P_hat = MX.zeros(Nx_p,Nx_p);
    
    % PREDICTION STEP
    z_hat = A_p*z_est(:,i-1) + B_p*u_p(:,i-1);
    P_hat = A_p*P(:,(i-2)*Nx_p+1:(i-1)*Nx_p)*A_p' + Q_p;
    z_est(:,i) = z_hat;

    % PREPROCESSING
    for m = 1:M
        % Get expected robot measurements with predicted states
        h_vec(m) = get_h(z_hat,x_rd(2*m-1,i),x_rd(2*m,i));

        % Compute Jacobian
        H(m,1) = h_vec(m) / z_hat(1);
        H(m,3) = h_vec(m) * (1/z_hat(3) - ((x_rd(2*m-1,i) - z_hat(5))^2 + (x_rd(2*m,i) - z_hat(7))^2));
        H(m,5) = h_vec(m) * (-2*z_hat(3)*(z_hat(5) - x_rd(2*m-1,i)));
        H(m,7) = h_vec(m) * (-2*z_hat(3)*(z_hat(7) - x_rd(2*m,i)));
    end
    
    % UPDATE STEP
    K = (P_hat*H')*inv(H*P_hat*H' + R_p);
    P(:,(i-1)*Nx_p+1:i*Nx_p) = P_hat - K*H*P_hat;
    
    % Define p as the diagonal elements of P
    p(:,i) = diag(P(:,(i-1)*Nx_p+1:i*Nx_p));

    % Collection :)
    z_hat_dum(:,i) = z_hat;
    P_hat_dum(:,(i-1)*Nx_p+1:i*Nx_p) = P_hat;
    K_dum(:,(i-1)*M+1:i*M) = K;
    H_dum(:, (i-1)*Nx_p+1:i*Nx_p) = H;

end

% Now that p has been defined in terms of the robot positions, we can write
% up the cost
cost_KF = 0;
for i = 2:Hp+1
    %cost_KF = cost_KF + p(:,i)'*Q*p(:,i);
    cost_KF = cost_KF + p(:,i)'*Q_vec;
end



power_cons = @(v) 0.0001*v^2 + 0.05;

% Constraints
barrier_dist = MX.zeros(M,Hp+1);
powers = MX.zeros(M,Hp+1);
robot_dist = MX.zeros(M,Hp+1);
e = MX.zeros(M,Hp+1);              % Energy of robot
e(:,1) = e0;
cost_slack = 0;
for i = 2:Hp+1
    % Robot dynamics
    if i <= Hu+1
        opti.subject_to(x_rd(:,i) == A_rd*x_rd(:,i-1) + B_rd*du(:,i-1));
    else
        opti.subject_to(x_rd(:,i) == A_rd*x_rd(:,i-1)); % DU(k > Hu) = 0
    end

    % Robot constraints
    for m = 1:M
        % Velocity and box constraints
        x_vel = x_rd(2*M+2*m-1,i);
        y_vel = x_rd(2*M+2*m,i);
        opti.subject_to(x_vel^2 + y_vel^2 <= vmax^2);

        % Box constraints
        opti.subject_to(xmin <= x_rd(2*m-1,i));
        opti.subject_to(xmax >= x_rd(2*m-1,i));
        opti.subject_to(ymin <= x_rd(2*m,i));
        opti.subject_to(ymax >= x_rd(2*m,i));

        % % Energy dynamics
        % x_pos = x_rd(2*m-1,i-1);
        % y_pos = x_rd(2*m,i-1);
        % powers(m,i-1) = power(x_pos,y_pos,sqrt(x_vel^2 + y_vel^2));
        % e(m,i) = e(m,i-1) + powers(m,i-1)*Ts;
        % 
        % % General energy constraints
        % opti.subject_to(0 <= e(m,i) + e_low_slack(m,i));   % Energy must not go under 0
        % opti.subject_to(e(m,i) <= 1 + e_high_slack(m,i));   % Energy must not exceed 1
        % 
        % % Energy constaints (barrier)
        % if(i == Hp+1)
        %     x_pos = x_rd(2*m-1,i);
        %     y_pos = x_rd(2*m,i);
        %     robot_dist(:,i) = sqrt((x_pos - charger_x)^2 + (y_pos - charger_y)^2);
        %     barrier_dist(m,i) = vmax*(e(m,i)/power_cons(vmax)); % Distance that robot m can get by going full speed in one direction
        %     opti.subject_to((x_pos - charger_x)^2 + (y_pos - charger_y)^2 <= barrier_dist(m,i)^2 + b_slack(m,i));
        % end
        % cost_slack = cost_slack + b_slack(m,i) + e_low_slack(m,i) + e_high_slack(m,i);
        
        % Slack constraints
        opti.subject_to(b_slack(m,i) >= 0);
        opti.subject_to(e_low_slack(m,i) >= 0);
        opti.subject_to(e_high_slack(m,i) >= 0);
    end
end

opti.set_initial(b_slack, ones(M,Hp+1));
opti.set_initial(e_low_slack, ones(M,Hp+1));
opti.set_initial(e_high_slack, ones(M,Hp+1));

% Define MPC 'object'
opti.minimize(cost + cost_KF/1000 + 1000*cost_slack);
solver_opts = struct();
%solver_opts.ipopt.max_iter = 5000;
%solver_opts.ipopt.tol = 1e-3;
solver_opts.ipopt.print_level = verbose_opt;
solver_opts.print_time = false;
opti.solver('ipopt', solver_opts);

%% Call MPC

% Warm start
if do_warm_start
    opti.set_initial(sol_prev.value_variables());
end


%sol = opti.solve();

% Solve
try
    sol = opti.solve();
    sol.value(x_rd);
    sol.value(e)
    %sol.value(p)
    % sol.value(opti.f())
    % sol.stats.iter_count
    fprintf("Objective: %.2f\n", sol.value(opti.f()))
    fprintf("Iterations: %d\n", sol.stats.iter_count)
catch
    disp("MPC FAILED!!! :((((")
    debug = opti.debug
    disp("")
end

% Outputs
X_opt = sol.value(x_rd(1:2*M,:));
U_opt = sol.value(x_rd(2*M+1:end,:));
P_trace = 0;


end