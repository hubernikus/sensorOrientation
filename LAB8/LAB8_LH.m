%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Sensor Orientation - EPFL
%           LAB 8
%
%           Author: Huber Lukas
%           Date: 2017-11-25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all; clear variables; clc;

%addpath(genpath('../functions'))
addpath(genpath('functions'))

%% Rescaling sensor error
clc;
g = 9.81; % [m/s^2] - Gravity

%% Ex 1
% Simulation Parameters
samplingFreq = 1; % [Hz]
dt_gps = 1/samplingFreq;

% Kalman filter frequency
dt_kf_list = [1,0.1];

% Initial conditions
omega0 = pi/100;     % Angular Rate [rad/s]
r_circ = 25;       % [m]

azim_init = 90/180*pi; % [rad]
phi0 = 90/180*pi; % [rad]
x0 = [0;r_circ];
vel0 = [-omega0*r_circ;0]; % initial velocity
initHeading = [phi0+azim_init];

%% Measurement Simulation
close all; 
[acc, gyro, x_real, v_real, phi_real, time] = ... 
               IMUsens_simulation(samplingFreq, omega0, r_circ, phi0, azim_init);

% GPS Simulation 
simga_x_GPS = 0.5; simga_y_GPS = 0.5; % [m]
x_simu = x_real + [whiteNoiseGen(size(x_real,2),1,simga_x_GPS); whiteNoiseGen(size(x_real,2),1,simga_y_GPS)];

%%
constVel = false;

for iter = 1:2
    dt_kf = dt_kf_list(iter);
    
    if constVel
        dim = 4;

        % A priori Info Conditions
        sigma_x = 10;% [m]
        sigma_v = 0.10;% [m]

        F = [0 0 1 0; 
             0 0 0 1;
             0 0 0 0;
             0 0 0 0];

        G = [0 0;
             0 0;
             1 0;
             0 1];

        Phi = eye(dim) + dt_kf*F;

        sigma_vt = 0.035; % Uncertainty in Model
        a = 0; % uniform constant motion
        sigma_at = 0;
        qvn = sigma_vt^2; qve = sigma_vt^2;
        W = diag([qvn, qve]);


        R = diag([simga_x_GPS^2, simga_y_GPS^2]);


        % State 
        %x = [r_circ; 0; x_simu(1,1);x_simu(2,1)];
        % x = [x, y, v_x, v_y];
        x_filt = [0; r_circ; -omega0*r_circ; 0];
        
        P0 = diag([sigma_x^2, sigma_x^2, sigma_v^2, sigma_v^2]);

        H = [1 0 0 0;
             0 1 0 0];
    else
        % Kalman Filtering
        dim = 6;

        % A priori Info Conditions
        sigma_p = 10;       % [m]
        sigma_v = 0.10;     % [m]
        sigma_a = 0.10;     % [m]

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

        Phi0 = eye(6) ...
            + [zeros(4,2), eye(4)*dt_kf; zeros(2,6)] ...
            + [zeros(2,4), eye(2)*dt_kf^2/2; zeros(4,6)];

        % Uncertainty of motion model
        sigma_at = 0.001; 
        a_dot = 0; % uniform constant acceleration
        qvn = sigma_at^2; qve = sigma_at^2;
        W = diag([qvn, qve]);

        % R [l x l] - Covariance matrix of measurment
        R = diag([simga_x_GPS^2, simga_y_GPS^2]);

        % x [n x 1]- State vector
        % x = [x, y, v_x, v_y, a_x, a_y];
        x_filt = [0; r_circ; -omega0*r_circ; 0; 0; -omega0^2*r_circ];

        % z [l x 1] - Measurements 

        % H [l x n]- Design Matrix ( z = H*x + v )
        H = [1 0 0 0 0 0;
             0 1 0 0 0 0];

        % P [n x n]- Covariance Matrix of the state vector
        P0 = diag([sigma_p^2, sigma_p^2, sigma_v^2, sigma_v^2, sigma_a^2, sigma_a^2]);
    end
    P = P0;

    % Initialization
    % Auxiliary matrix A
    A = [-F, G*W*G'; ...
         zeros(dim), F'] * dt_kf;

    B = expm(A);

    Phi = B(dim+1:2*dim,dim+1:2*dim)';
    Q_k = Phi*B(1:dim,dim+1:2*dim);

    % Prediction
    x_tild(:,1) = Phi*x_filt(:,1);
    x_tilde = Phi*x_filt(:,1);
    P_tilde = Phi*P*Phi' + Q_k;

    N_sim = size(x_simu,2);
    
    %N_sim = 30
    for ii = 1:N_sim
        % Mesurement --- z, R ----
        
        % Gain (weight)
        K = P_tilde*H'/(H*P_tilde*H' + R);

        % State update
        innovation(:,iter,ii) = (x_simu(:,ii)- H*x_tild(:,1+((ii-1)*dt_gps/dt_kf) ));
        x_filt(:,ii+1) = x_tild(:,1+((ii-1)*dt_gps/dt_kf)) + K*innovation(:,iter,ii);
     
        % Covariance update
        P = (eye(dim)- K*H)*P_tilde;
        
        % KF-predicted
        sigma_pred(iter,ii) = sqrt(sum(diag(P)));

        % Prediction - Inner Loop
        % Output ---- x_est, P ---

        % Auxiliary matrix A
        A = [-F, G*W*G'; ...
             zeros(dim), F'] * dt_gps;

        B = expm(A);

        Phi = B(dim+1:2*dim,dim+1:2*dim)';
        Q_k = Phi*B(1:dim,dim+1:2*dim);
        
        x_tilde = Phi*x_filt(:,ii+1);
        %x_tilde = 0;
        
        for jj = 1:dt_gps/dt_kf    
            
            % Output ---- x_est, P ---
            
            % Auxiliary matrix A
            A = [-F, G*W*G'; ...
                 zeros(dim), F'] * dt_kf;

            B = expm(A);

            Phi = B(dim+1:2*dim,dim+1:2*dim)';
            Q_k = Phi*B(1:dim,dim+1:2*dim);

            % Prediction
            if(jj==1) % first round take measurement
                x_tild(:,jj+1 + (ii-1)*dt_gps/dt_kf) = Phi*x_filt(:,ii+1);
                P_tilde = Phi*P*Phi' + Q_k;
            else    % take prediction after
                x_tild(:,jj+1 + (ii-1)*dt_gps/dt_kf) = Phi*x_tild(:,jj+(ii-1)*dt_gps/dt_kf);
                P_tilde = Phi*P_tilde*Phi' + Q_k;
            end
        end
    end
%     
    % Empirical standard deviaion characterizing GPS
    sigma_gps_x = std(x_simu(1,:)-x_real(1,:));
    sigma_gps_y = std(x_simu(2,:)-x_real(2,:));
    sigma_gps(iter) = sqrt(sigma_gps_x^2+sigma_gps_y^2)

    % Empirical standard deviaion characterizing GPS
    sigma_filt_x = std(x_filt(1,1:end-1)-x_real(1,:));
    sigma_filt_y = std(x_filt(2,1:end-1)-x_real(2,:));
    sigma_filt(iter) = sqrt(sigma_filt_x^2+sigma_filt_y^2)
    
    % Empirical standard deviaion characterizing GPS
    meanSqr_gps(iter) = sum(sigma_gps_x^2+sigma_gps_y^2)

    % Empirical standard deviaion characterizing GPS
    meanSqr_filt(iter) = sum(sigma_filt_x^2+sigma_filt_y^2)
    
%     figure;
%     subplot(2,1,1)
%     plot(x_simu(1,:)-x_real(1,:),'r'); hold on;
%     plot(x_filt(1,1:end-1)-x_real(1,:),'g'); 
%     
%     subplot(2,1,2)
%     plot(x_simu(2,:)-x_real(2,:),'r'); hold on;
%     plot(x_filt(2,1:end-1)-x_real(2,:),'g');
   x_KF(:,:,iter) =  x_filt;
end


%%
% close all;

% Rotational plot
figure('Position',[100,200,800,700]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.0)
plot(x_real(1,1:N_sim),x_real(2,1:N_sim),'b-x'); hold on;
plot(x_simu(1,1:N_sim),x_simu(2,1:N_sim),'ro')
%plot(x_simu(1,1:4),x_simu(2,1:4),'bo')
iter = 1;
plot(x_KF(1,:, iter),x_KF(2,:, iter),'g-x')
iter = 2;
plot(x_KF(1,:, iter),x_KF(2,:, iter),'m-x')
ylabel('North [m]'); xlabel('East [m]')
legend('Real measurement', 'Measurement','Kalman Filter 1 Hz', 'Kalman Filter 0.1 Hz')

print(sprintf('fig/kalmanFilter_sigma%1.4f.eps',sigma_at),'-depsc')


%% I. Improvement in the  Positioning accuracy & Anticipate

%sabilisationTime = sabilisationTime
sigma_pred_fin = sigma_pred(:,end)

% Print Latex Table
fileID = fopen('table_STD.tex','w');
fprintf(fileID,'Real GPS Positioning [m]& %3.6f & %3.6f  \\\\ \\hline \n',sigma_gps );
fprintf(fileID,'Filtered Positioning [m]& %3.6f & %3.6f  \\\\ \\hline \n',sigma_filt );
fclose(fileID);

%% II. Anticipated accuracy (compare to const model)
%
figure('Position',[100,200,800,400]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
ylabel('Predicted positioning')
for iter = 1:size(sigma_pred,1)
    plot(dt_gps*(0:size(sigma_pred,2)-1),sigma_pred(iter,:)); hold on;
end
xlabel('Time [s]')
legend('Prediction interval: 1 Hz','Prediction interval: 0.1 Hz')

if constVel
    print('fig/KFpredict_constVel','-depsc')
else
    print('fig/KFpredict','-depsc')
end


%% III. Velocity distribution error
time = [0:size(v_real,2)-1]*dt_gps;

figure('Position',[100,200,800,700]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
subplot(2,1,1)
plot(time,x_KF(3,1:end-2,1),'r'); hold on;
plot(time,x_KF(3,1:end-2,2),'g'); hold on;
plot(time,v_real(1,:),'b'); hold on;
ylabel('Velocity East [m/s]')
subplot(2,1,2)
plot(time,x_KF(4,1:end-2,1),'r'); hold on;
plot(time,x_KF(4,1:end-2,2),'g'); hold on;
plot(time,v_real(2,:),'b'); hold on;
ylabel('Velocity North [m/s]')
xlabel('Time [s]')
print('fig/distribution_velocity','-depsc')

figure('Position',[100,200,800,700]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
subplot(2,1,1)
plot(time,x_KF(3,1:end-2,1)-v_real(1,:),'r'); hold on;
plot(time,x_KF(3,1:end-2,2)-v_real(1,:),'g'); hold on;
ylabel('Velocity East [m/s]')
subplot(2,1,2)
plot(time,x_KF(4,1:end-2,1)-v_real(2,:),'r'); hold on;
plot(time,x_KF(4,1:end-2,2)-v_real(2,:),'g'); hold on;
ylabel('Velocity North [m/s]')
xlabel('Time [s]')
if constVel
    print('fig/distribution_velocityError_constVel','-depsc')
else
    print('fig/distribution_velocityError','-depsc')
end
%% IV. Histogramm innovation  - 
time = [0:size(v_real,2)-1]*dt_gps;

figure('Position',[100,200,800,700]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
subplot(2,1,1)
ii = 1; dimen = 1;
jj = 2;
hist([squeeze(innovation(dimen,ii,1:end-2)), ...
      squeeze(innovation(dimen,jj,1:end-2))]); hold on;
std_hist = std([squeeze(innovation(dimen,ii,1:end-2)); ...
      squeeze(innovation(dimen,jj,1:end-2))]);
mean_hist = mean([squeeze(innovation(dimen,ii,1:end-2)); ...
      squeeze(innovation(dimen,jj,1:end-2))]);
N_vals = 100;
histMax = 50;
xLims = xlim;
xVals = linspace(xLims(1),xLims(2), N_vals);
f_norm = histMax*1/(std_hist*sqrt(2*pi)).*exp(-0.5*(xVals-mean_hist).^2/std_hist);
plot(xVals, f_norm, 'k--')

xlabel('Innovation East [m/s]')
legend('Prediction interval: 1 Hz','Prediction interval: 0.1 Hz')
subplot(2,1,2)
ii = 1; dimen = 2;
jj = 2;
hist([squeeze(innovation(dimen,ii,1:end-2)), ...
      squeeze(innovation(dimen,jj,1:end-2))]); hold on;
std_hist = std([squeeze(innovation(dimen,ii,1:end-2)); ...
      squeeze(innovation(dimen,jj,1:end-2))]);
mean_hist = mean([squeeze(innovation(dimen,ii,1:end-2)); ...
      squeeze(innovation(dimen,jj,1:end-2))]);
N_vals = 100;
histMax = 50;
xLims = xlim;
xVals = linspace(xLims(1),xLims(2), N_vals);
f_norm = histMax*1/(std_hist*sqrt(2*pi)).*exp(-0.5*(xVals-mean_hist).^2/std_hist);
plot(xVals, f_norm, 'k--')
xlabel('Innovation North [m/s]')

if sigma_at == 0.001
    print('fig/histogramm_innovation_all_sigma_at0001','-depsc')
elseif sigma_at == 0.1
    print('fig/histogramm_innovation_all_sigma_at01','-depsc')
else
    print('fig/histogramm_innovation_all','-depsc')
end

%% IV. Histogramm Innovation 
stabTime = 30;
figure('Position',[100,200,800,700]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
subplot(2,1,1)
ii = 1; dimen = 1;
jj = 2;
hist([squeeze(innovation(dimen,ii,stabTime:end)), ...
      squeeze(innovation(dimen,jj,stabTime:end))]); hold on;
std_hist = std([squeeze(innovation(dimen,ii,stabTime:end)); ...
      squeeze(innovation(dimen,jj,stabTime:end))]);
mean_hist = mean([squeeze(innovation(dimen,ii,stabTime:end)); ...
      squeeze(innovation(dimen,jj,stabTime:end))]);
N_vals = 100;
histMax = 50;
xLims = xlim;
xVals = linspace(xLims(1),xLims(2), N_vals);
f_norm = histMax*1/(std_hist*sqrt(2*pi)).*exp(-0.5*(xVals-mean_hist).^2/std_hist);
plot(xVals, f_norm, 'k--')

xlabel('Innovation East [m/s]')
legend('Prediction interval: 1 Hz','Prediction interval: 0.1 Hz')
subplot(2,1,2)
ii = 1; dimen = 2;
jj = 2;
hist([squeeze(innovation(dimen,ii,stabTime:end)), ...
      squeeze(innovation(dimen,jj,stabTime:end))]); hold on;
std_hist = std([squeeze(innovation(dimen,ii,stabTime:end)); ...
      squeeze(innovation(dimen,jj,stabTime:end))]);
mean_hist = mean([squeeze(innovation(dimen,ii,stabTime:end)); ...
      squeeze(innovation(dimen,jj,stabTime:end))])
N_vals = 100;
histMax = 50;
xLims = xlim;
xVals = linspace(xLims(1),xLims(2), N_vals);
f_norm = histMax*1/(std_hist*sqrt(2*pi)).*exp(-0.5*(xVals-mean_hist).^2/std_hist);
plot(xVals, f_norm, 'k--')
xlabel('Innovation North [m/s]')

if sigma_at == 0.001
    print('fig/histogramm_innovation_cut_sigma_at0001','-depsc')
elseif sigma_at == 0.1
    print('fig/histogramm_innovation_cut_sigma_at01','-depsc')
else
    print('fig/histogramm_innovation_cut','-depsc')
end


%% V. Changes position 
if not(constVel)
    time = [0:size(x_filt,2)-1]*dt_gps;
    figure('Position',[100,200,900,900]);
    set(groot,'DefaultAxesFontSize',14)
    set(groot,'DefaultLineLineWidth',1.2)
    subplot(3,1,1)
    iter = 1;
    plot(time,x_KF(1,:,iter),'r'); hold on;
    plot(time,x_KF(2,:,iter),'r--'); hold on;
    iter = 2;
    plot(time,x_KF(1,:,iter),'g'); hold on;
    plot(time,x_KF(2,:,iter),'g--'); hold on;
    ylabel('Position [m]')

    subplot(3,1,2)
    iter = 1;
    plot(time,x_KF(3,:,iter),'r'); hold on;
    plot(time,x_KF(4,:,iter),'r--'); hold on;
    iter = 2;
    plot(time,x_KF(3,:,iter),'g'); hold on;
    plot(time,x_KF(4,:,iter),'g--'); hold on;
    ylabel('Velocity [m/s]')
    legend('Up. freq. 1 Hz (East)','Up. freq. 1 Hz (North)', ...
        'Up. freq. 0.1 Hz (East)','Up. freq. 0.1 Hz (North)')

    subplot(3,1,3)
    iter = 1;
    plot(time,x_KF(5,:,iter),'r'); hold on;
    plot(time,x_KF(6,:,iter),'r--'); hold on;
    iter = 2;
    plot(time,x_KF(5,:,iter),'g'); hold on;
    plot(time,x_KF(6,:,iter),'g--'); hold on;
    ylabel('Acceleration [m/s^2]')
    xlabel('Time [s]')
    if sigma_at == 0.001
        print('fig/kalmanFilter_subplots_vel_acc_sigma_at0001','-depsc')
    elseif sigma_at == 0.1
        print('fig/kalmanFilter_subplots_vel_acc_sigma_at01','-depsc')
    else
        print('kalmanFilter_subplots_vel_acc_','-depsc')
    end
    
    time = [0:size(x_tild,2)-1]*dt_kf;
    figure('Position',[100,200,800,700]);
    set(groot,'DefaultAxesFontSize',14)
    set(groot,'DefaultLineLineWidth',1.2)
    subplot(3,1,1)
    plot(time,x_tild(1,:),'r'); hold on;
    plot(time,x_tild(2,:),'g'); hold on;
    ylabel('Position [m]')

    subplot(3,1,2)
    plot(time,x_tild(3,:),'r'); hold on;
    plot(time,x_tild(4,:),'g'); hold on;
    ylabel('Velocity [m/s]')
    legend('East','North')

    subplot(3,1,3)
    plot(time,x_tild(5,:),'r'); hold on;
    plot(time,x_tild(6,:),'g'); hold on;
    ylabel('Acceleration [m/s^2]')
    xlabel('Time [s]')
    print('fig/kalmanFilter_tilde_subplots','-depsc')
end

%% VI. Decrease increase velocity



%%
figure('Position',[100,200,800,600]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
subplot(2,1,1)
for iter = 1:size(innovation,2)
    plot(dt_gps*(0:size(innovation,3)-1),squeeze(innovation(1,iter,:))); hold on;
end
ylabel('Innovation (east) [m]')
legend('Prediction interval: 1 Hz','Prediction interval: 0.1 Hz')
subplot(2,1,2)
for iter = 1:size(innovation,2)
    plot(dt_gps*(0:size(innovation,3)-1),squeeze(innovation(1,iter,:))); hold on;
end
ylabel('Innovation (north) [m]')
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
ylabel('Innovation (east) [m]')
legend('Prediction interval: 1 Hz','Prediction interval: 0.1 Hz')
subplot(2,1,2)
for iter = 1:size(innovation,2)
    plot(dt_gps*(sabilisationTime(iter)-1:size(innovation,3)-1),squeeze(innovation(2,iter,sabilisationTime(iter):end))); hold on;
end
ylabel('Innovation (north) [m]')
xlabel('Time [s]')
print('fig/stableInovations','-depsc')

