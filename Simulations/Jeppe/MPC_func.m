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
du = opti.variable(M*Nu_r, Hu+1);   % Delta u (change in robot control input)  
x = opti.variable(M*Nx_r + M*Nu_r,Hp+1);            % Robot positions defined as optimization variables but are bounded by constraints  

% States and varialbes (NOT to be optimized over)
z_est = opti.parameter(Nx_p,Hp+1);          % Kalman filter states
P = opti.parameter(Nx_p,(Hp+1)*Nx_p);       % Kalman filter error covariance matrix
p = opti.parameter(Nx_p,Hp+1);              % KF error covariance diagonal elements
z_hat = opti.parameter(Nx_p,Hp+1);
P_hat = opti.parameter(Nx_p,(Hp+1)*Nx_p);
%H = opti.parameter(M,(Hp+1)*Nx_p);
%h_vec = opti.parameter(M,Hp+1);
%K = opti.parameter(Nx_p,(Hp+1)*M);

% Enforce initial conditions
opti.subject_to(du(:,1) == zeros(M*Nu_r,1))
opti.subject_to(x(1:M*Nx_r,1) == x0);
opti.set_value(z_est(:,1), z_est_0);
opti.set_value(P(:,1:Nx_p), P_0);
opti.set_value(p(:,1), zeros(Nx_p,1));  % This is to ensure that indices are matched
opti.set_value(z_hat(:,1), zeros(Nx_p,1));% This is to ensure that indices are matched
opti.set_value(P_hat(:,1:Nx_p), zeros(Nx_p,Nx_p));
%opti.set_value(H(:,1:Nx_p), zeros(M,Nx_p));
%opti.set_value(h_vec(:,1), zeros(M,1));
%opti.set_value(K(:,1:M), zeros(Nx_p,M));

% % Decision variables (to be optimized over)
% DU = opti.variable(M*l_a,Hu);
% Za = opti.variable(M*n_a,Hp+1);
% 
% P_KF_est = opti.variable(8*Hp,8);
% P_KF_pred = opti.variable(8*(Hp+1),8);
% 
% Z_KF = opti.variable(8,Hp+1);
% h_KF = opti.variable(M,Hp);
% 
% H_KF = opti.variable(M,8);
% K_KF = opti.variable(8,M);
% 
% % Parameters (NOT to be optimized over)
% Z0 = opti.parameter(M*n,1);
% Uprev = opti.parameter(M*l,1);
% Z0_KF = opti.parameter(8,1);
% P0_KF = opti.parameter(8,8);


%% Create Cost Function
% Control input
cost = 0;
for i = 2:Hu+1
    cost = cost + du(:,i)'*R*du(:,i);
end

% Kalman filter iterations
% z_est = MX.zeros(Nx_p,1);
% P = MX.zeros(Nx_p,Nx_p);
for i = 2:Hp+1
    % Define temp. variables
    H = MX.zeros(M,Nx_p);
    h_vec = MX.zeros(M,1);
    K = MX.zeros(Nx_p,M);
    z_hat = MX.zeros(Nx_p,1);
    P_hat = MX.zeros(Nx_p,Nx_p);
    
    % PREDICTION STEP
    z_hat = A_p*z_est(:,i-1) + B_p*u_p(:,i-1);
    P_hat = A_p*P(:,(i-2)*Nx_p+1:(i-1)*Nx_p)*A_p' + Q_p;

    % PREPROCESSING
    for m = 1:M
        % Get expected robot measurements with predicted states
        h_vec(m) = get_h(z_hat,x(2*m-1,i),x(2*m));

        % Compute Jacobian
        H(m,1) = h_vec(m) / z_est(1,i);
        H(m,3) = h_vec(m) * (1/z_est(3,i) - ((x(2*m-1,i) - z_est(5,i))^2 + (x(2*m) - z_est(7,i))^2));
        H(m,5) = h_vec(m) * (-2*z_est(3,i)*(z_est(5,i) - x(2*m-1,i)));
        H(m,7) = h_vec(m) * (-2*z_est(3,i)*(z_est(7,i) - x(2*m)));
    end
    
    % UPDATE STEP
    K = (P_hat*H') / (H*P_hat*H' + R_p);
    P(:,(i-1)*Nx_p+1:i*Nx_p) = P_hat - K*H*P_hat;
    
    % Define p as the diagonal elements of P
    p(:,i) = diag(P(:,(i-1)*Nx_p+1:i*Nx_p));
end

% Now that p has been defined in terms of the robot positions, we can write
% up the cost
cost_KF = 0;
for i = 1:Hp
    cost_KF = cost_KF + p(:,i)'*Q*p(:,i);
end

% Constraints
for i = 2:Hp+1
    % Robot dynamics
    if i <= Hu+1
        opti.subject_to(x(:,i) == A_aN*Za(:,k) + B_aN*DU(:,k));
    else
        opti.subject_to(Za(:,k+1) == A_aN*Za(:,k)); % DU(k > Hu) = 0
    end
end

% Constraints
Za0 = [reshape(Z0,n,M); reshape(Uprev,l,M)];
opti.subject_to(Za(:,1) == reshape(Za0,M*(n+l),1));
for k = 1:Hp
    % Dynamics
    if k <= Hu
        opti.subject_to(Za(:,k+1) == A_aN*Za(:,k) + B_aN*DU(:,k));
    else
        opti.subject_to(Za(:,k+1) == A_aN*Za(:,k)); % DU(k > Hu) = 0
    end

    % Box constraints
    opti.subject_to(xmin <= Za(1:n_a:M*n_a,k+1) <= xmax);
    opti.subject_to(ymin <= Za(2:n_a:M*n_a,k+1) <= ymax);

    % Velocity constraints
    for j = 1:M
        x_vel = Za(3+(j-1)*n_a,k+1);
        y_vel = Za(4+(j-1)*n_a,k+1);
        opti.subject_to(x_vel^2 + y_vel^2 <= vmax^2);
    end

    % % Collision avoidance
    % for i = 1:M
    %     for j = i+1:M
    %         xi = Za((i-1)*n_a + 1, k+1);
    %         yi = Za((i-1)*n_a + 2, k+1);
    %         xj = Za((j-1)*n_a + 1, k+1);
    %         yj = Za((j-1)*n_a + 2, k+1);
    % 
    %         dist_squared = (xi - xj)^2 + (yi - yj)^2;
    %         opti.subject_to(dist_squared >= dmin^2);
    %     end
    % end
end

% Outputs
Z_opt = [];
U_opt = [];
for i = 1:n_a:M*n_a
    Z_opt = [Z_opt; Za(i,:); Za(i+1,:)];
    U_opt = [U_opt; Za(i+2,:); Za(i+3,:)];
end

% Define MPC 'object'
opti.minimize(lambda1*cost + lambda2*cost_KF);

solver_opts = struct;
solver_opts.ipopt.max_iter = 2000;
solver_opts.ipopt.tol = 1e-6;
solver_opts.ipopt.print_level = verbose_opt;
% solver_opts.print_time = false;
opti.solver('ipopt', solver_opts);

%% Call MPC
% Initialization
Z0_val = rob_init(:,1); % Initial robot states
Uprev_val = rob_init(:,2); % Previous control input
%Z0_KF_val = kf_init(:,1); % Initial prediction for kalman filter
%P0_KF_val = kf_init(:,2:end); % Initial error covariance
Z0_KF_val = z_est(:,1);
P0_KF_val = P(:,:,1);

opti.set_value(Z0, Z0_val);
opti.set_value(Uprev, Uprev_val);
opti.set_value(Z0_KF, Z0_KF_val);
opti.set_value(P0_KF, P0_KF_val);

% Solve
sol = opti.solve();

% Outputs
Z_sol = sol.value(Z_opt);
U_sol = sol.value(U_opt);

P_sol = sol.value(P_KF_pred);
P_trace = [];
for k = 1:Hp+1
    P_trace = [P_trace, log(abs(trace(P_sol(8*k-7:8*k,:))))];
end

end

