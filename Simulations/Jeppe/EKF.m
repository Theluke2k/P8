function [x,P] = EKF(x_prev,P_prev,A_prev,B_prev,u_prev,Q_prev, y, R, h, H, pos_rob)
    % INPUTS
    % x_prev  : States from k-1
    % P_prev  : Error covariance matrix from k-1
    % A_prev  : System matrix for process dynamics from k-1
    % B_prev  : Input matrix for process dynamics from k-1
    % u_prev  : Input vector from from k-1
    % Q_prev  : State noise covariance matrix from k-1
    % y       : Measurements from k-1
    % R       : Measurement noise covariance matrix from k-1
    % h       : Function to compute the measurements (h vector in EKF)
    % H       : Function to comp. Jacobian of h with states and robot pos.
    % pos_rob : State vector for robots (robot positions)
    % 
    % OUTPUTS
    % x  : State vector for time k
    % P  : Error covariance matrix for time k

    % PREDICTION STEP
    x_hat = A_prev*x_prev + B_prev*u_prev;
    P_hat = A_prev*P_prev*A_prev' + Q_prev;

    % UPDATE STEP
    v = y - h(x_hat, pos_rob);
    K = P_hat*H(x_hat,pos_rob)'*inv(H(x_hat,pos_rob)*P_hat*H(x_hat,pos_rob)' + R);
    x = x_hat + K*v;
    P = P_hat - K*H(x_hat,pos_rob)*P_hat;
end

