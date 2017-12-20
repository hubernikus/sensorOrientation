function [x_filt, innovation, sigma_pred] = kalmanFilter(x_simu, P0, x0, H, Phi, Q, R, F,G)

x_filt = x0;

P = P0;

dim = size(x0,1);

% initialize variables
sigma_pred = zeros(1,size(x_simu,2)+1);
innovation = zeros(dim,size(x_simu,2)+1);

for i = 1:size(x_simu,2)-1
    % Prediction
    x_tild = Phi*x_filt(:,i);
    P_tilde = Phi*P*Phi' + Q;

    % Mesurement --- z, R ----

    % Gain (weight)
    K = P_tilde*H'/(H*P_tilde*H' + R);

    % State update
    innovation(:,i) = (x_s  imu(:,i+1)- H*x_tild);
    x_filt(:,i+1) = x_tild + K*innovation(:,i);

    % Covariance update
    %P = (eye(dim)- K*H)*P_tilde;
    P = (eye(dim)- K*H)*P_tilde;

    % KF-predicted
    sigma_pred(i) = sqrt(sum(diag(P)));

    % Output ---- x_est, P ---
end

end