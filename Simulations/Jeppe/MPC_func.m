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
h_p = kf_init{8};          % Function handle to compute vector h
H_p = kf_init{9};          % Function handle to compute Jacobian matrix H
error_cv_selec = kf_init{10};% Selection of diagonal elements to have in cost function

% Useful numbers
Nx_p = size(A_p,1);       % Number of states in process model
Nu_p = size(B_p,2);       % Number of inputs in process model

% Define initial conditions
z_est = zeros(Nx_p,Hp+1);      % Estimated process states
P = zeros(Nx_p,Nx_p,Hp+1);        % KF error covariance matrix
z_est(:,1) = kf_init{1};    % KF state from k-1
P(:,:,1) = kf_init{2};      % KF error covariance matrix fom k-1

%% Robot model and Tuning
% Tuning for a single robot
Q_single = [1 0;
     0 1];
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
Q = kron(eye(Nx_p), Q_single);
R = kron(eye(Nu_p), R_single);

%% MPC object
% Setup opti
opti = Opti();

% Decision variables (to be optimized over)
du = opti.variable(M*Nu, Hu);   % Delta u (change in robot control input)

% States and varialbes (NOT to be optimized over)
p = opti.variable(Nx_p,Hp+1);  % KF error covariance diagonal elements
x = opti.variable(M*Nx,Hp+1);   % Robot positions defined as optimization variables but are bounded by constraints    

% Enforce initial conditions on variables
opti.subject_to(x(:,1) == rob_init{1});

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
for k = 1:Hp
    cost = cost + du(:,k)'*Q*du(:,k); % Assuming R = I
end

% Kalman Filter
dummy1 = zeros(Nx_p);
for k = 1:Hp
    
end


cost_KF = 0;
Z_KF(:,1) = Z0_KF;
P_KF_pred(1:8,:) = P0_KF; % P(k|k-1)
for k = 1:Hp
    for j = 1:M
        % Positions
        x_pos = Za(1+(j-1)*n_a,k);
        y_pos = Za(2+(j-1)*n_a,k);

        % Measurement predictions
        h_KF(j,k) = h(Z_KF(:,k), x_pos, y_pos);

        % Linearize
        H_KF(j,1) = h_KF(j,k)/(Z_KF(1,k)+eps);
        H_KF(j,3) = -h_KF(j,k)*((x_pos-Z_KF(5,k))^2 + (y_pos-Z_KF(7,k))^2) + h_KF(j,k)/(Z_KF(3,k)+eps);
        H_KF(j,5) = -2*h_KF(j,k)*Z_KF(3,k)*(Z_KF(5,k)-x_pos);
        H_KF(j,7) = -2*h_KF(j,k)*Z_KF(3,k)*(Z_KF(7,k)-y_pos);
    end
    K_KF = P_KF_pred(8*k-7:8*k,:)*H_KF'*inv(H_KF*P_KF_pred(8*k-7:8*k,:)*H_KF' + R_KF);

    P_KF_est(8*k-7:8*k,:) = (eye(8) - K_KF*H_KF)*P_KF_pred(8*k-7:8*k,:)*(eye(8) - K_KF*H_KF)' + K_KF*R_KF*K_KF'; % P(k|k)

    % Time update
    Z_KF(:,k+1) = A_KF*Z_KF(:,k) + B_KF*u_KF(:,k);
    P_KF_pred(8*(k+1)-7:8*(k+1),:) = A_KF*P_KF_est(8*k-7:8*k,:)*A_KF' + Q;

    % Add to cost
    cost_KF = cost_KF + P_KF_est(8*k-7,1)/20  + P_KF_est(8*k-5,3)*100 + P_KF_est(8*k-3,5) + P_KF_est(8*k-1,7);
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

