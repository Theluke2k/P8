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
N = n_robots; % Number of robots

%% KALMAN FILTER

% KF Dynamic model
tau_beta = 0.998;       %spread increase parameter (exponential decay of beta)

% Inputs
beta0 = 0; %steady-state value for beta
u_beta = beta0*(1-tau_beta)/dt; % constant input to achieve beta0 in steady-state
u_M = 2; % Leakage rate offset
u_KF = ones(2,Hp);
u_KF(1,:) = u_KF(1,:)*u_M;
u_KF(2,:) = u_KF(2,:)*u_beta;

A_KF = [1 dt 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0;
     0 0 tau_beta dt 0 0 0 0;
     0 0 0 1 0 0 0 0;
     0 0 0 0 1 dt 0 0;
     0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 dt;
     0 0 0 0 0 0 0 1];

B_KF = zeros(8,2);
B_KF(1,1) = dt; 
B_KF(3,2) = dt;

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
R_KF = sigma_measurement*eye(N); % measurement noise covariance

%% Robot model

% Robot model (integrator)
A = eye(2); % z = [x;y]
B = Ts*eye(2); % u = [vx;vy]
C = eye(2);
D = 0;
n = size(A,1);
l = size(B,2);
m = size(C,1);

% Reformulation to Delta u
A_a = [A,B; zeros(n,n),eye(n)]; % X(k) = [x(k),y(k),vx(k-1),vy(k-1)]^T
B_a = [B; eye(l)]; % Du(k) = [vx(k)-vx(k-1), vy(k)-vy(k-1)]^T
C_a = [C, zeros(n,m)]; % za(k) = [x(k),y(k),vx(k-1),vy(k-1)]^T
D_a = 0; % DU(k) = [vx(k)-vx(k-1), vy(k)-vy(k-1)]^T
n_a = size(A_a,1);
l_a = size(B_a,2);
m_a = size(C_a,1);

% Expanded to include all N robots (DEFINED INCORRECTLY -> Assumes [x1;y1;vx1;vy1;x2;y2;...] should assume [x1;y1;x2;y2;...;vx1;vy1;...])
A_aN = kron(eye(N), A_a);
B_aN = kron(eye(N), B_a);
C_aN = kron(eye(N), C_a);
D_aN = 0;

%% MPC object
% Setup opti
opti = Opti();

% Decision variables (to be optimized over)
DU = opti.variable(N*l_a,Hu);
Za = opti.variable(N*n_a,Hp+1);

P_KF_est = opti.variable(8*Hp,8);
P_KF_pred = opti.variable(8*(Hp+1),8);

Z_KF = opti.variable(8,Hp+1);
h_KF = opti.variable(N,Hp);

H_KF = opti.variable(N,8);
K_KF = opti.variable(8,N);

% Parameters (NOT to be optimized over)
Z0 = opti.parameter(N*n,1);
Uprev = opti.parameter(N*l,1);
Z0_KF = opti.parameter(8,1);
P0_KF = opti.parameter(8,8);

% Cost function
cost = 0;
for k = 1:Hu
    cost = cost + DU(:,k)'*DU(:,k); % Assuming R = I
end

% KF algorithm
cost_KF = 0;
Z_KF(:,1) = Z0_KF;
P_KF_pred(1:8,:) = P0_KF; % P(k|k-1)
for k = 1:Hp
    for j = 1:N
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
Za0 = [reshape(Z0,n,N); reshape(Uprev,l,N)];
opti.subject_to(Za(:,1) == reshape(Za0,N*(n+l),1));
for k = 1:Hp
    % Dynamics
    if k <= Hu
        opti.subject_to(Za(:,k+1) == A_aN*Za(:,k) + B_aN*DU(:,k));
    else
        opti.subject_to(Za(:,k+1) == A_aN*Za(:,k)); % DU(k > Hu) = 0
    end

    % Box constraints
    opti.subject_to(xmin <= Za(1:n_a:N*n_a,k+1) <= xmax);
    opti.subject_to(ymin <= Za(2:n_a:N*n_a,k+1) <= ymax);

    % Velocity constraints
    for j = 1:N
        x_vel = Za(3+(j-1)*n_a,k+1);
        y_vel = Za(4+(j-1)*n_a,k+1);
        opti.subject_to(x_vel^2 + y_vel^2 <= vmax^2);
    end

    % % Collision avoidance
    % for i = 1:N
    %     for j = i+1:N
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
for i = 1:n_a:N*n_a
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
Z0_KF_val = kf_init(:,1); % Initial prediction for kalman filter
P0_KF_val = kf_init(:,2:end); % Initial error covariance

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

