function [X_filt, innovation, sigma_pred, X_tild] ...
            = kalmanFilter_differential(model, x_init, dx_init, x_gps, meas_imu, time_gps, time_imu, err)
% Extract variables from model
P = model.P0;
F = model.F;
G = model.G;
dimX = model.dimZ;
dimZ = model.dimX;
R = model.R;
W = model.W;
H_k = model.H;
Phi = model.Phi0;
Q_k = model.Q0;

% Initialize matrices
N_gps = size(time_gps,2);
N_imu = size(time_imu,2);

dX_filt = zeros(dimZ, N_gps) ;
X_filt = zeros(dimZ, N_gps) ;
innovation = zeros(dimX,N_gps);
sigma_pred = zeros(dimZ, N_gps);
dX_tild = zeros(dimZ,N_imu);
X_tild = zeros(dimZ,N_imu);

% First values 
X_filt(:,1) = x_init;
X_tild(:,1) = x_init;

dX_filt(:,1) = dx_init;
dX_tild(:,1) = dx_init;

% Initialization
% Prediction
P_tilde = Phi*P*Phi' + Q_k;

sigma_pred(1) = sqrt(sum(diag(P)));

% Iteration counter variables
it_gps = 1;
it_imu = 2;
firstInnerLoop = true;

figure('Position',[100 100 800 800])
%plot_gps = plot(x_gps(2,1:it_gps),x_gps(1,1:it_gps),'bo'); 
plot_gps = plot(x_gps(2,:),x_gps(1,:),'bo'); 
hold on; grid on; axis equal;
%   xlim([-600,600]); ylim([-600,600]);
plot_filt = plot(X_filt(5,1:it_gps),X_filt(4,1:it_gps),'g-x');
plot_tild = plot(X_tild(5,1:it_imu),X_tild(4,1:it_imu),'r-x');
legend('GPS simulation','Filtered estimation','Prediction with IMU')

%while(time_gpsn(it_gps) <= timeFinal)
while(it_gps <= N_gps)
    % Mesurement --- z, R ----

    % Gain (weight)
    K_k = P_tilde*H_k'/(H_k*P_tilde*H_k' + R);

    % State update
    innovation(:,it_gps) = (x_gps(:,it_gps)- H_k*X_tild(:,it_imu));
    
    dX_filt(:, it_gps+1) = dX_tild(:,it_imu) + K_k*innovation(:,it_gps);
    X_filt(:, it_gps+1) = X_filt(:, it_gps) + dX_filt(:, it_gps+1);
    X_tild(:, it_imu) = X_filt(:, it_gps) + dX_filt(:, it_gps+1); % Replace X_tild, with filtered state
    
    set(plot_filt,'XData',X_filt(5,1:it_gps+1),'YData',X_filt(4,1:it_gps+1));

    % Covariance update
    P = (eye(dimZ)- K_k*H_k)*P_tilde;

    % KF-predicted
    sigma_pred(it_gps+1) = sqrt(sum(diag(P)));

    % Prediction - Inner Loop
    % Output ---- x_est, P ---
    alpha = X_filt(1,it_gps+1);

    firstInnerLoop = true; % first Time doing the inner loop
    while(time_imu(it_imu) <= time_gps(it_gps) ) % bigger or bigger equal???
        % Output ---- x_est, P ---
        
        % Prediction
        if(firstInnerLoop) % first round take GPS comparison

        end
        R_3 = [cos(alpha), -sin(alpha); sin(alpha),cos(alpha)];
        acc_m = R_3*meas_imu(2:3,it_imu)
        
        F11 = [0 0 0 0 0; 
               -acc_m(1) 0 0 0 0;
               acc_m(2) 0 0 0 0;
               0 1 0 0 0;
               0 0 1 0 0];

        F21 = zeros(4,5);

        F12 = [1 1 0 0; 
               0 0 cos(alpha) -sin(alpha);
               0 0 sin(alpha) cos(alpha);
               0 0 0 0;
               0 0 0 0];

        F22 = diag([0, -err.gyro_GM1.beta, -err.acc1_GM1.beta, -err.acc2_GM1.beta]);

        model.F = [F11,F12; F21,F22];

        G11 = [1 0 0;
               0 cos(alpha) -sin(alpha);
               0 sin(alpha) cos(alpha);
               0 0 0;
               0 0 0];

        G22 = [0 0 0;
               1 0 0;
               0 1 0;
               0 0 1];

        model.G = blkdiag(G11,G22);
        
        
        % Auxiliary matrix A
        A = [-F, G*W*G'; ...
             zeros(dimZ), F'] * (time_imu(it_imu+1)-time_imu(it_imu));

        B = expm(A);

        Phi = B(dimZ+1:2*dimZ,dimZ+1:2*dimZ)';
        Q_k = Phi*B(1:dimZ,dimZ+1:2*dimZ);

        % Prediction
        if(firstInnerLoop) % first round take measurement
            dX_tild(:,it_imu+1) = Phi*dX_filt(:,it_gps+1);
            P_tilde = Phi*P*Phi' + Q_k;

        else    % take prediction after
            dX_tild(:,it_imu+1) = Phi*dX_tild(:,it_imu);
            P_tilde = Phi*P_tilde*Phi' + Q_k;
            
        end
        X_tild(:, it_imu+1) = X_tild(:, it_imu) + dX_tild(:, it_imu+1);
                
        alpha = X_tild(1,it_imu);
        firstInnerLoop = false;
        it_imu = it_imu+1;
        set(plot_tild,'XData',X_tild(5,1:it_imu),'YData',X_tild(4,1:it_imu));
    end
   %set(plot_gps,'XData',x_gps(2,1:it_gps-1),'YData',x_gps(1,1:it_gps-1));
    % increment counter
   it_gps = it_gps+1;

end
end
