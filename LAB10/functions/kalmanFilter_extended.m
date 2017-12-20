function [x_filt, innovation, sigma_pred, x_tild] ...
            = kalmanFilter_extended(model, x_init, x_gps, meas_imu, time_gps, time_imu, err)
% Extract variables from model
P = model.P0;
F = model.F;
G = model.G;
dim = model.dimPhi;
R = model.R;
W = model.W;
H = model.H;
Phi = model.Phi;
Q_k = model.Q_k;

% Initialize matrices
N_gps = size(time_gps,2);
N_imu = size(time_imu,2);

x_filt = zeros(dim, N_gps) ;
innovation = zeros(dim,N_gps);
sigma_pred = zeros(dim, N_gps);
x_tild = zeros(dim,N_imu);

% Initialization
% Prediction
x_tild(:,1) = Phi*x_filt(:,1);
P_tilde = Phi*P*Phi' + Q_k;

sigma_pred(1) = sqrt(sum(diag(P)));

% Assign first value
x(:,1) = x_init;

% Iteration counter variables
it_gps = 1;
it_imu = 1;

% Get final time
timeFinal = time_gps(end);

while(time_gps(it_gps) <= timeFinal)
    % Mesurement --- z, R ----

    % Gain (weight)
    K = P_tilde*H'/(H*P_tilde*H' + R);

    % State update
    innovation(:,it_gps) = (x_gps(:,it_gps)- H*x_tild(:,it_imu));
    
    x_filt(:,it_gps+1) = x_tild(:,it_imu) + K*innovation(:,it_gps);
    
    % Covariance update
    P = (eye(dim)- K*H)*P_tilde;

    % KF-predicted
    sigma_pred(it_gps+1) = sqrt(sum(diag(P)));

    % Prediction - Inner Loop
    % Output ---- x_est, P ---

    % Auxiliary matrix A

    while(time_imu(it_imu) < time_gps(it_gps) ) % bigger or bigger equal???
        % Output ---- x_est, P ---
        
        % Prediction
        if(jj==1) % first round take measurement
            alpha = x_filt(1,it_gps+1);
        else    % take prediction after
            alpha = x_tild(1,it_imu);
        end


        F11 = [0 0 0 0 0; 
               0 0 0 0 0;
               0 0 0 0 0;
               0 0 0 0 1;
               0 0 0 1 0];

        F21 = zeros(4,5);

        F12 = [1 1 0 0; 
               0 0 cos(alpha) sin(alpha);
               0 0 sin(alpha) cos(alpha);
               0 0 0 0;
               0 0 0 0];

        F22 = diag([0, -err.gyro_GM1.beta, -err.acc1_GM1.beta, -err.acc2_GM1.beta]);

        model.F = [F11,F12; F21,F22];

        G11 = [1 0 0;
               0 cos(alpha) sin(alpha);
               0 0 0;
               0 0 0];

        G22 = [0 0 0;
               0 0 0;
               1 0 0;
               0 1 0;
               0 0 1];

        model.G = blkdiag(G11,G22);
        
        
        % Auxiliary matrix A
        A = [-F, G*W*G'; ...
             zeros(dim), F'] * dt_kf;

        B = expm(A);

        Phi = B(dim+1:2*dim,dim+1:2*dim)';
        Q_k = Phi*B(1:dim,dim+1:2*dim);

        % Prediction
        if(jj==1) % first round take measurement
            x_tild(:,it_imu+1) = Phi*x_filt(:,it_gps+1);
            P_tilde = Phi*P*Phi' + Q_k;
        else    % take prediction after
            x_tild(:,it_imu+1) = Phi*x_tild(:,it_imu);
            P_tilde = Phi*P_tilde*Phi' + Q_k;
        end
        
        % increment counter
        it_imu = it_imu+1;
    end
    % increment counter
    it_gps = it_gps+1;
end
    
end
