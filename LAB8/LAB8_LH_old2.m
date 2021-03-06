%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Sensor Orientation - EPFL
%           LAB 7 
%
%           Author: Huber Lukas
%           Date: 2017-11-18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all; clear variables; clc;

%addpath(genpath('../functions'))
addpath(genpath('functions'))

%% Rescaling sensor errorc
clc;
g = 9.81; % [m/s^2] - Gravity

%% Ex 1
% Simulation Parameters
samplingFreq = 1/2; % [Hz]
dt_gps = 1/samplingFreq;

% Simulation Parameters
samplingFreq = 1; % [Hz]
%dt_gps = 1/samplingFreq;

% Kalman filter frequency
dt_kf = 1;
dt_kf = 0.1;

% Initial conditions
omega0 = pi/100;     % Angular Rate [rad/s]
r_circ = 25;       % [m]

azim_init = 90/180*pi; % [rad]
phi0 = 90/180*pi; % [rad]
x0 = [0;r_circ];
vel0 = [-omega0*r_circ;0]; % initial velocity
initHeading = [phi0+azim_init];

%% Measurement Simulation
for iter = 1:1

    [acc, gyro, x_real, v_real, phi_real, time] = ... 
               IMUsens_simulation(samplingFreq, omega0, r_circ, phi0, azim_init);

    % GPS Simulation 
    simga_x_GPS = 0.5; simga_y_GPS = 0.5; % [m]
    x_simu = x_real + [whiteNoiseGen(size(x_real,2),1,simga_x_GPS); whiteNoiseGen(size(x_real,2),1,simga_y_GPS)];

    % Kalman Filtering
    close all; 

    dim = 6;

    % A priori Info Conditions
    sigma_p = 10;% [m]
    sigma_v = 0.10;% [m]
     sigma_a = 0.10;% [m]
     
    F = [0 0 1 0 0 0; 
         0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         0 0 0 0 0 0;
         0 0 0 0 0 0];

    G = [0 0;
         0 0;
         0 0;
         0 0;
         1 0;
         0 1];


    Phi = eye(6) ...
        + [zeros(4,2), eye(4)*dt_kf; zeros(2,6)] ...
        + [zeros(2,4), eye(2)*dt_kf^2/2; zeros(4,6)];

    sigma_at = 0.01; % Uncertainty of motion model
    a_dot = 0; % uniform constant acceleration
    qvn = sigma_at^2; qve = sigma_at^2;
    W = diag([qvn, qve]);

    %Q = diag([sigma_x^2, sigma_x^2, sigma_v^2, sigma_v^2]);
%         Q = [1/3*qvn*dt_gps^3, 0 , 1/2*dt_gps^2*qvn, 0;
%          0, 1/3*qve*dt_gps^3, 0 , 1/2*dt_gps^2*qve;
%          1/2*qvn*dt_gps^2, 0 , dt_gps*qvn, 0 
%          0, 1/2*qve*dt_gps^2, 0 , dt_gps*qve];
    % R [l x l] - Covariance matrix of measurment
    R = diag([simga_x_GPS^2, simga_y_GPS^2]);
    
    % x [n x 1]- State vector
    % x = [x, y, v_x, v_y, a_x, a_y];
    x_filt = [0; r_circ; -omega0*r_circ; 0; 0; 0];

    % z [l x 1] - Measurements 
        
    % H [l x n]- Design Matrix ( z = H*x + v )
    H = [1 0 0 0 0 0;
         0 1 0 0 0 0];
    
    % P [n x n]- Covariance Matrix of the state vector
    P0 = diag([sigma_p^2, sigma_p^2, sigma_v^2, sigma_v^2, sigma_a^2, sigma_a^2]);
    P = P0;
    
    % Initialization
    % Auxiliary matrix A
    A = [-F, G*W*G'; ...
         zeros(dim), F'] * dt_kf;

    B = expm(A);

    Phi = B(dim+1:2*dim,dim+1:2*dim)';
    Q = Phi*B(1:dim,dim+1:2*dim);
    
    % Prediction
    x_tild = Phi*x_filt(:,1);
    P_tilde = Phi*P*Phi' + Q;

    for i = 1:size(x_simu,2)-1

        % Mesurement --- z, R ----

        % Gain (weight)
        K = P_tilde*H'/(H*P_tilde*H' + R);

        % State update
        innovation(:,iter,i) = (x_simu(:,i+1)- H*x_tild);
        x_filt(:,i+1) = x_tild + K*innovation(:,iter,i);

        % Covariance update
        %P = (eye(dim)- K*H)*P_tilde;
        P = (eye(dim)- K*H)*P_tilde;

        % KF-predicted
        sigma_pred(iter,i) = sqrt(sum(diag(P)));
        
        % Prediction inner loop
        for jj = 1:dt_gps/dt_kf
            % Output ---- x_est, P ---
            % Auxiliary matrix A
            A = [-F, G*W*G'; ...
                 zeros(dim), F'] * dt_kf;

            B = expm(A);

            Phi = B(dim+1:2*dim,dim+1:2*dim)';
            Q = Phi*B(1:dim,dim+1:2*dim);
            % Prediction
            
            if(jj==1)
                x_tild = Phi*x_filt(:,i+1);
            else
                x_tild = Phi*x_tild;
            end
            P_tilde = Phi*P*Phi' + Q;
        end


    end
    
    % Empirical standard deviaion characterizing GPS
    sigma_gps_x = std(x_simu(1,:)-x_real(1,:))
    sigma_gps_y = std(x_simu(2,:)-x_real(2,:))
    sigma_gps(iter) = sqrt(sigma_gps_x^2+sigma_gps_y^2)

    % Empirical standard deviaion characterizing GPS
    sigma_filt_x = std(x_filt(1,:)-x_real(1,:))
    sigma_filt_y = std(x_filt(2,:)-x_real(2,:))
    sigma_filt(iter) = sqrt(sigma_filt_x^2+sigma_filt_y^2);
end

%%
close all;

figure('Position',[100,200,800,700]);
plot(x_real(1,:),x_real(2,:),'b'); hold on;
plot(x_simu(1,:),x_simu(2,:),'ro')
%plot(x_simu(1,1:4),x_simu(2,1:4),'bo')
plot(x_filt(1,:),x_filt(2,:),'g-x')
ylabel('North [m]'); xlabel('East [m]')
legend('Real measurement', 'Measurement','Kalman Filter')
print('fig/kalmanFilter','-depsc')


figure('Position',[100,200,800,400]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
ylabel('Predicted positioning')
for iter = 1:size(sigma_pred,1)
    plot(dt_gps*(0:size(sigma_pred,2)-1),sigma_pred(iter,:)); hold on;
end
xlabel('Time [s]')
print('fig/KFpredict','-depsc')



figure('Position',[100,200,800,600]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
subplot(2,1,1)
ylabel('Innovation (east)')
for iter = 1:size(innovation,2)
    plot(dt_gps*(0:size(innovation,3)-1),squeeze(innovation(1,iter,:))); hold on;
end
subplot(2,1,2)
for iter = 1:size(innovation,2)
    plot(dt_gps*(0:size(innovation,3)-1),squeeze(innovation(1,iter,:))); hold on;
end
ylabel('Innovation (north)')
xlabel('Time [s]')
print('fig/allInovations','-depsc')

%%
figure('Position',[100,200,800,600]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
subplot(2,1,1)
for iter = 1:size(innovation,2)
    sabilisationTime(iter) = find(abs((sigma_pred(iter,(2:end))-sigma_pred(iter,(1:end-1)))...
                        ./sigma_pred(iter,(2:end))) < 1e-4 ,1);
    plot(dt_gps*(sabilisationTime(iter)-1:size(innovation,3)-1),squeeze(innovation(1,iter,sabilisationTime(iter):end))); hold on;
end
ylabel('Innovation (east)')
subplot(2,1,2)
for iter = 1:size(innovation,2)
    plot(dt_gps*(sabilisationTime(iter)-1:size(innovation,3)-1),squeeze(innovation(2,iter,sabilisationTime(iter):end))); hold on;
end
ylabel('Innovation (north)')
xlabel('Time [s]')
print('fig/stableInovations','-depsc')


sabilisationTime = sabilisationTime
sigma_pred_fin = sigma_pred(:,end)
%% Print Latex Table


fileID = fopen('table_STD.tex','w');
fprintf(fileID,'Real GPS Positioning [m]& %3.4f & %3.4f & %3.4f & %3.4f & %3.4f & %3.4f  \\\\ \\hline \n',sigma_gps, mean(sigma_gps) );
fprintf(fileID,'Filtered Positioning [m]& %3.4f & %3.4f & %3.4f & %3.4f & %3.4f & %3.4f  \\\\ \\hline \n',sigma_filt, mean(sigma_filt) );
fclose(fileID);