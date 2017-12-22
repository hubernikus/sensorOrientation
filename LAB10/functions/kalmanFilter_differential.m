function [X_filt, innovation, sigma_pred, X_tild] ...
            = kalmanFilter_differential(model, x_init, dx_init, x_gps, meas_imu, time_gps, time_imu, err)
% Extract variables from model
P = model.P0;
F = model.F;
G = model.G;
dimX = model.dimX;
dimZ = model.dimZ;
R = model.R;
W = model.W;
H_k = model.H;
Phi = model.Phi0;
Q_k = model.Q0;

% Initialize matrices
N_gps = size(time_gps,2);
N_imu = size(time_imu,2);

gyro = meas_imu(1,:);
acc = meas_imu(2:3,:);


dX_filt = zeros(dimX, N_gps) ;
X_filt = zeros(dimX, N_gps) ;

innovation = zeros(dimZ,N_gps);
sigma_pred = zeros(dimX, N_gps);

dX_tild = zeros(dimX,N_imu);
X_tild = zeros(dimX,N_imu);
X_star = zeros(dimX,N_imu);


% First values 
X_filt(:,1) = x_init;
X_tild(:,1) = x_init;
X_star(:,1) = x_init;

dX_filt(:,1) = dx_init;
dX_tild(:,1) = dx_init;

% Initialization
% Prediction
P_tilde = Phi*P*Phi' + Q_k;

sigma_pred(1) = sqrt(sum(diag(P)));

% Iteration counter variables
it_gps = 1;
it_imu = 1;
firstInnerLoop = true;

figure('Position',[100 100 800 800])
%plot_gps = plot(x_gps(2,1:it_gps),x_gps(1,1:it_gps),'bo'); 
plot_gps = plot(x_gps(2,:),x_gps(1,:),'bo'); 
hold on; grid on; axis equal;
ylim([400,600]); xlim([-50,150]);
plot_filt = plot(X_filt(5,1:it_gps),X_filt(4,1:it_gps),'g-x');
plot_tild = plot(X_tild(5,1:it_imu),X_tild(4,1:it_imu),'r-x');
legend('GPS simulation','Filtered estimation','Prediction with IMU')

%while(time_gpsn(it_gps) <= timeFinal)
while(it_gps <= N_gps-1)
    % Mesurement --- z, R ----

    % Gain (weight)
    K_k = P_tilde*H_k'/(H_k*P_tilde*H_k' + R);

    % State update
    dZ = (x_gps(:,it_gps)- H_k*X_tild(:,it_imu));
    innovation(:,it_gps) = (dZ- H_k*dX_tild(:,it_imu));
    %innovation(:,it_gps) = (dZ- 0);
    
    dX_filt(:, it_gps+1) = dX_tild(:,it_imu) + K_k*innovation(:,it_gps);
    %dX_filt(:, it_gps+1) = 0 + K_k*innovation(:,it_gps);
    X_filt(:, it_gps+1) = X_tild(:, it_imu) + dX_filt(:, it_gps+1) ...
                                *(time_imu(it_imu+1)-time_imu(it_imu));
    %dX_filt(:, it_gps+1)  = 0;
    
    %X_star(:, it_imu) = X_star(:, it_gps) + dX_filt(:, it_gps+1); % Replace X_tild, with filtered state
    
    set(plot_filt,'XData',X_filt(5,1:it_gps+1),'YData',X_filt(4,1:it_gps+1));
    X_tild(:, it_imu) = X_filt(:, it_gps+1); % Replace X_tild, with filtered state
    
    % Covariance update
    P = (eye(dimX)- K_k*H_k)*P_tilde;

    % KF-predicted
    sigma_pred(it_gps+1) = sqrt(sum(diag(P)));

    % Prediction - Inner Loop
    % Output ---- x_est, P ---
    
    firstInnerLoop = true; % first Time doing the inner loop
    while(time_imu(it_imu) <= time_gps(it_gps+1) ) % bigger or bigger equal???
        % Output ---- x_est, P ---
        
        % Prediction x_dot
        X_star(:,it_imu+1) = X_star(:,it_imu);
        [X_star(1:5,it_imu+1)] = ...
            inertialNavigationStep(X_star(1:5,it_imu), ...
                                    acc(:,it_imu), gyro(it_imu), (time_imu(it_imu+1)-time_imu(it_imu)));
        
        alpha = X_star(1,it_imu+1);
        R_3 = [cos(alpha), -sin(alpha); sin(alpha),cos(alpha)];
        acc_m = R_3*acc(:,it_imu);
        
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
             zeros(dimX), F'] * (time_imu(it_imu+1)-time_imu(it_imu));

        B = expm(A);

        Phi = B(dimX+1:2*dimX,dimX+1:2*dimX)';
        Q_k = Phi*B(1:dimX,dimX+1:2*dimX);

        % Prediction
        if(firstInnerLoop) % first round take measurement
            dX_tild(:,it_imu+1) = Phi*dX_filt(:,it_gps+1);
            P_tilde = Phi*P*Phi' + Q_k;

        else    % take prediction after
            dX_tild(:,it_imu+1) = Phi*dX_tild(:,it_imu);
            P_tilde = Phi*P_tilde*Phi' + Q_k;
            
        end
        X_tild(:, it_imu+1) = X_filt(:, it_gps) + dX_tild(:, it_imu+1);
                
        it_imu = it_imu+1;
        set(plot_tild,'XData',X_tild(5,1:it_imu),'YData',X_tild(4,1:it_imu));
        firstInnerLoop = false;
    end
   %set(plot_gps,'XData',x_gps(2,1:it_gps-1),'YData',x_gps(1,1:it_gps-1));
    % increment counter
   it_gps = it_gps+1;

end
end
