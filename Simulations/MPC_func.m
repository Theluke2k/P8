function [Z_sol,U_sol,P_trace,opti] = MPC_func(rob_init, kf_init, mpc_params, cost_params, n_robots, map_bounds, min_dist, sim_params, verbose_opt)
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

% Define initial conditions
%z_est = zeros(Nx_p,Hp+1);      % Estimated process states
%P = zeros(Nx_p,Nx_p,Hp+1);        % KF error covariance matrix


%% Robot model and Tuning
% Tuning for a single robot
Q = eye(Nx_p);
R_single = [1 0;
     0 1];

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

x0 = rob_init{1};

%% MPC object
% Setup opti
opti = Opti();

% Decision variables (to be optimized over)
du = opti.variable(M*Nu_r, Hu);   % Delta u (change in robot control input)  

% States and varialbes (NOT to be optimized over)
x = opti.parameter(M*Nx_r,Hp+1);            % Robot positions defined as optimization variables but are bounded by constraints  
z_est = opti.parameter(Nx_p,Hp+1);          % Kalman filter states
P = opti.parameter(Nx_p,(Hp+1)*Nx_p);       % Kalman filter error covariance matrix
p = opti.parameter(Nx_p,Hp+1);              % KF error covariance diagonal elements
z_hat = opti.parameter(Nx_p,Hp+1);
P_hat = opti.parameter(Nx_p,(Hp+1)*Nx_p);

% Enforce initial conditions on variables
opti.set_value(x(:,1), x0);
opti.set_value(z_est(:,1), z_est_0);
opti.set_value(P(:,1:Nx_p), P_0);
opti.set_value(p(:,1), zeros(Nx_p,1));  % This is to ensure that indices are matched
opti.set_value(z_hat(:,1), zeros(Nx_p,1));% This is to ensure that indices are matched
opti.set_value(P_hat(:,1:Nx_p), zeros(Nx_p,Nx_p));



%% Create Cost Function
% Control input
cost = 0;
for i = 1:Hu
    cost = cost + du(:,i)'*R*du(:,i);
end

% Error covariance trace
cost_KF = 0;
for i = 1:Hp
    cost_KF = cost_KF + p(:,i)'*Q*p(:,i);
end

% Kalman filter iterations
H = MX.zeros(M,(Hp+1)*Nx_p);
h_vec = MX.zeros(M,Hp+1);
K = MX.zeros(Nx_p,(Hp+1)*M);
for i = 2:Hp+1
    Pcols = (i-1)*Nx_p+1:i*Nx_p;
    Hcol = (i-1)*Nx_p+1;
    Hcols = Hcol:Hcol+Nx_p-1;
    Kcols = (i-1)*M+1:i*M;

    % PREDICTION STEP
    z_hat(:,i) = A_p*z_est(:,i-1) + B_p*u_p(:,i-1)
    P_hat(:,Pcols) = A_p*P(:,Pcols-Nx_p)*A_p' + Q_p;
    
    % PREPROCESSING
    for m = 1:M
        h_vec(m,i) = get_h(z_hat(:,i),x(2*m-1,i),x(2*m));
        
        % Compute Jacobian
        H(m,Hcol) = h_vec(m,i) / z_est(1,i);
        H(m,Hcol+2) = h_vec(m,i) * (1/z_est(3,i) - ((x(2*m-1,i) - z_est(5,i))^2 + (x(2*m) - z_est(7,i))^2));
        H(m,Hcol+4) = h_vec(m,i) * (-2*z_est(3,i)*(z_est(5,i) - x(2*m-1,i)));
        H(m,Hcol+6) = h_vec(m,i) * (-2*z_est(3,i)*(z_est(7,i) - x(2*m)));
    end
    
    % UPDATE STEP
    K(:,Kcols) = (P_hat(:,Pcols)*H(:,Hcols)') / (H(:,Hcols)*P_hat(:,Pcols)*H(:,Hcols)' + R_p);
    P(:,Pcols) = P_hat(:,Pcols) - K(:,Kcols)*H(:,Hcols)*P_hat(:,Pcols);
end

end

