%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Sensor Orientation - EPFL
%           LAB 10
%
%           Author: Huber Lukas
%           Date: 2017-11-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('functions'))

clear all; clc
close all;

% Change values to update all figures
figSave = false;
figName ='';

%% Rescaling sensor error
g = 9.81; % [m/s^2] - Gravity

%% Simulation Parameters
freq_IMU = 100; % [Hz]
dt_kf = 1/freq_IMU;

%samplingFreq = 1; % [Hz]

% Kalman filter frequency
t_start = 0;
%dt_kf_list = [1,0.1];

% Initial conditions
omega0 = pi/100;     % Angular Rate [rad/s]
r_circ = 500;       % [m]

azim_init = 90/180*pi; % [rad]
phi0 = 90/180*pi; % [rad]
initHeading = [phi0+azim_init];


%% 1 - Generate reference trajectory 
%  2 - Simulate 2D - IMU measurement
pos_0 = [r_circ,0];
vel_0 = [0, omega0*r_circ];

[acc, gyro, x_real, v_real, phi_real, time_imu] = ...
            IMUsens_simulation(freq_IMU, omega0, r_circ, phi0, azim_init);

%% 2 - Add noise to IMU measurement
[acc_simu, gyro_simu, init_errorVariables] ...
                    = IMUsens_errors(acc, gyro, freq_IMU);

figure('Position',[50,50,800,1000]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
subplot(3,1,1)
plot(time_imu,gyro_simu, 'b'); hold on; grid on;
ylabel('Gyroscope [rad/s]')
xlim([time_imu(1),time_imu(end)])
subplot(3,1,2)
plot(time_imu,acc_simu(1,:), 'b'); hold on; grid on; 
ylabel('Acceleromter north [m/s^2]')
xlim([time_imu(1),time_imu(end)])
subplot(3,1,3)
plot(time_imu,acc_simu(2,:), 'b'); hold on; grid on;
xlim([time_imu(1),time_imu(end)])
ylabel('Acceleromter east[m/s^2]')
xlabel('Time [s]')
grid on;
close all;

if figSave; print(strcat('fig/acc_gyro_simu',figName),'-depsc'); end

%% 3 - Simualte GPS signale 
dt_gps = 2;                          % [s]
freq_GPS = 1/dt_gps;

simga_x_gps = 1.0; simga_y_gps = 1.0; % [m]
sigma_gps = [simga_x_gps,simga_y_gps];

% GPS Simulation 
indGPS = 1:freq_IMU/freq_GPS:size(x_real,2);

pos_gps = x_real(:,indGPS); 
time_gps = time_imu(:,indGPS); 

pos_gps_simu = pos_gps + [whiteNoiseGen(size(pos_gps,2),1,simga_x_gps); 
                                whiteNoiseGen(size(pos_gps,2),1,simga_y_gps)];

%% 4 - Kalman Initialization
p_GPS_0 = pos_gps_simu(:,1);

[model, X_init, dX_init] ...
            = init_kalmanFilter(dt_kf, r_circ, omega0,[gyro_simu(1);acc_simu(:,1)], sigma_gps, p_GPS_0, init_errorVariables);

%% 4 - Kalman filter loop

[x_filt, innovation, sigma_pred, x_tild] ...
            = kalmanFilter_differential(model, X_init, dX_init, pos_gps_simu, [gyro_simu;acc_simu], time_gps, time_imu, init_errorVariables);
%%
% Rotational plot
figure('Position',[100,200,800,700]);
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.0)
plot(pos_gps_simu(2,:),pos_gps_simu(1,:),'bo'); hold on;
plot(x_filt(5,:),x_filt(4,:),'g-x'); hold on;
plot(x_tild(5,:),x_tild(4,:),'r-x');
ylabel('North [m]'); xlabel('East [m]')
legend('Real measurement','Kalman Filter', 'Prediction')
axis equal;
grid on;

if figSave; print(strcat('fig/kalmanFilter_sigma',figName),'-depsc'); end

%%
% 
% %% Kalman Filtering
% for ii = 1:5
% switch ii
%     case 1
%         % % Personalize figure name
%         fac_processingNoise = 1; % 0.1, 1, 10
%         modelType = 'circularMotion';
%         figName = modelType;
%     case 2
%         fac_processingNoise = 1; % 0.1, 1, 10
%         modelType = 'constAcc';
%         figName = modelType;
%     case 3
%         fac_processingNoise = 1; % 0.1, 1, 10
%         modelType = 'constVel';
%         figName = modelType;
%     case 4
%         fac_processingNoise = 100; % 0.1, 1, 10
%         modelType = 'circularMotion';
%         figName = '_procNoise10';
%     case 5
%         fac_processingNoise = 0.01; % 0.1, 1, 10
%         modelType = 'circularMotion';
%         figName = '_procNoise01';
% end
%  
% clear x_KF innovation sigma_pred x_simu
% 
% for iter = 1:2 % change to 2
%     %% Measurement Simulation
%     [acc, gyro, x_real, v_real, phi_real, time] = ... 
%                IMUsens_simulation(samplingFreq, omega0, r_circ, phi0, azim_init);
% 
%     % GPS Simulation 
%     simga_x_GPS = 0.5; simga_y_GPS = 0.5; % [m]
%     sigma_gps = [simga_x_GPS,simga_y_GPS];
%     x_simu(:,:,iter) = x_real + [whiteNoiseGen(size(x_real,2),1,simga_x_GPS); 
%                                 whiteNoiseGen(size(x_real,2),1,simga_y_GPS)];
% 
%     N_sim = size(x_simu,2);
% 
%     dt_kf = 0.1;
%     
%     % Kalman Filtering
%      [model, x_filt] = ...
%              init_kalmanFilter(modelType, dt_kf, r_circ, omega0, sigma_gps, fac_processingNoise);
%     
%     [x_KF(:,:,iter), innovation(:,:,iter), sigma_pred(:,iter), x_tild, x_filt_rad] ...
%                 = kalmanFilter_extended(model, x_simu(:,:,iter), x_filt, dt_kf, dt_gps);
% end
% 
% %%
% 
% % Rotational plot
% figure('Position',[100,200,800,700]);
% set(groot,'DefaultAxesFontSize',14)
% set(groot,'DefaultLineLineWidth',1.0)
% plot(x_real(2,1:N_sim),x_real(1,1:N_sim),'b-x'); hold on;
% iter = 1;
% plot(x_simu(2,1:N_sim,iter),x_simu(1,1:N_sim,iter),'ro'); hold on;
% iter = 2;
% plot(x_simu(2,1:N_sim,iter),x_simu(1,1:N_sim,iter),'go')
% iter = 1;
% plot(x_KF(2,:, iter),x_KF(1,:, iter),'r-x')
% iter = 2;
% plot(x_KF(2,:, iter),x_KF(1,:, iter),'g-x')
% ylabel('North [m]'); xlabel('East [m]')
% legend('Real measurement', 'Simulation 1', 'Simulation 2','Kalman Filter 1', 'Kalman Filter 2')
% axis equal;
% grid on;
% print(strcat('fig/kalmanFilter_sigma',figName),'-depsc')
% 
% %% I. Velocity distribution error
% time = [0:size(v_real,2)-1]*dt_gps;
% 
% figure('Position',[100,200,800,700]);
% set(groot,'DefaultAxesFontSize',14)
% set(groot,'DefaultLineLineWidth',1.2)
% subplot(2,1,1)
% plot(time,x_KF(3,1:end-2,1),'r'); hold on;
% plot(time,x_KF(3,1:end-2,2),'g'); hold on;
% plot(time,v_real(1,:),'b'); hold on;
% ylabel('Velocity North [m/s]')
% grid on;
% subplot(2,1,2)
% plot(time,x_KF(4,1:end-2,1),'r'); hold on;
% plot(time,x_KF(4,1:end-2,2),'g'); hold on;
% plot(time,v_real(2,:),'b'); hold on;
% ylabel('Velocity East [m/s]')
% xlabel('Time [s]')
% grid on;
% print('fig/distribution_velocity','-depsc')
% 
% figure('Position',[100,200,800,700]);
% set(groot,'DefaultAxesFontSize',14)
% set(groot,'DefaultLineLineWidth',1.2)
% subplot(2,1,1)
% plot(time,x_KF(3,1:end-2,1)-v_real(1,:),'r'); hold on;
% plot(time,x_KF(3,1:end-2,2)-v_real(1,:),'g'); hold on;
% ylabel('Velocity North [m/s]')
% grid on;
% subplot(2,1,2)
% plot(time,x_KF(4,1:end-2,1)-v_real(2,:),'r'); hold on;
% plot(time,x_KF(4,1:end-2,2)-v_real(2,:),'g'); hold on;
% ylabel('Velocity East [m/s]')
% xlabel('Time [s]')
% grid on;
% 
% print(strcat('fig/distribution_velocityError',figName),'-depsc')
% 
% 
% %% II. Improvement in the  Positioning accuracy & Anticipate
% stdGPS = empiricalSTD(x_simu, x_real)
% stdFilt = empiricalSTD(x_KF, x_real)
% 
% %sabilisationTime = sabilisationTime
% sigma_pred_fin = sigma_pred(:,end);
% 
% % Print Latex Table
% fileID = fopen(strcat('table_STD',figName,'.tex'),'w');
% fprintf(fileID,'Real GPS Positioning [m]& %3.6f %3.6f  \\\\ \\hline \n', stdGPS );
% fprintf(fileID,'Filtered Positioning [m]& %3.6f & %3.6f  \\\\ \\hline \n', stdFilt );
% fclose(fileID);
% 
% %% III. Anticipated accuracy (compare to const model)
% %
% figure('Position',[100,200,800,400]);
% set(groot,'DefaultAxesFontSize',14)
% set(groot,'DefaultLineLineWidth',1.2)
% ylabel('Predicted positioning')
% for iter = 1:size(sigma_pred,2)
%     plot(dt_gps*(0:size(sigma_pred,1)-1),sigma_pred(:,iter)); hold on;
% end
% xlabel('Time [s]')
% legend('Simulation 1','Simulation 2')
% grid on;
% 
% print(strcat('fig/KF_predictSTD',figName),'-depsc')
% 
% 
% %% IV. Histogramm innovation  - 
% time = [0:size(v_real,2)-1]*dt_gps;
% 
% figure('Position',[100,200,800,700]);
% set(groot,'DefaultAxesFontSize',14)
% set(groot,'DefaultLineLineWidth',1.2)
% subplot(2,1,1)
% ii = 1; dimen = 1;
% jj = 2;
% hist([squeeze(innovation(dimen,1:end-2,ii)); ...
%       squeeze(innovation(dimen,1:end-2,jj))]'); hold on;
% std_hist = std([squeeze(innovation(dimen,1:end-2,ii)), ...
%       squeeze(innovation(dimen,1:end-2,jj))]);
% mean_hist = mean([squeeze(innovation(dimen,1:end-2,ii)), ...
%       squeeze(innovation(dimen,1:end-2,jj))]);
% N_vals = 100;
% histMax = 50;
% xLims = xlim;
% xVals = linspace(xLims(1),xLims(2), N_vals);
% f_norm = histMax*1/(std_hist*sqrt(2*pi)).*exp(-0.5*(xVals-mean_hist).^2/std_hist);
% %plot(xVals, f_norm, 'k--')
% xlabel('Radius [m/s]')
% grid on;
% legend('Simulation 1','Simulation 2')
% 
% subplot(2,1,2)
% ii = 1; dimen = 2;
% jj = 2;
% hist([squeeze(innovation(dimen,1:end-2,ii)); ...
%       squeeze(innovation(dimen,1:end-2,jj))]'); hold on;
% std_hist = std([squeeze(innovation(dimen,1:end-2,ii)), ...
%       squeeze(innovation(dimen,1:end-2,jj))]);
% mean_hist = mean([squeeze(innovation(dimen,1:end-2,ii)), ...
%       squeeze(innovation(dimen,1:end-2,jj))]);
% N_vals = 100;
% histMax = 50;
% xLims = xlim;
% xVals = linspace(xLims(1),xLims(2), N_vals);
% f_norm = histMax*1/(std_hist*sqrt(2*pi)).*exp(-0.5*(xVals-mean_hist).^2/std_hist);
% %plot(xVals, f_norm, 'k--')
% xlabel('Orientation [m/s]')
% grid on;
% 
% print(strcat('fig/histogramm_innovation_all',figName),'-depsc')
% 
% %% IV. Histogramm Innovation 
% stabTime = 30;
% figure('Position',[100,200,800,700]);
% set(groot,'DefaultAxesFontSize',14)
% set(groot,'DefaultLineLineWidth',1.2)
% subplot(2,1,1)
% ii = 1; dimen = 1;
% jj = 2;
% hist([squeeze(innovation(dimen,1:end-2,ii)); ...
%       squeeze(innovation(dimen,1:end-2,jj))]'); hold on;
% std_hist = std([squeeze(innovation(dimen,1:end-2,ii)), ...
%       squeeze(innovation(dimen,1:end-2,jj))]);
% mean_hist = mean([squeeze(innovation(dimen,1:end-2,ii)), ...
%       squeeze(innovation(dimen,1:end-2,jj))]);
% N_vals = 100;
% histMax = 50;
% xLims = xlim;
% xVals = linspace(xLims(1),xLims(2), N_vals);
% f_norm = histMax*1/(std_hist*sqrt(2*pi)).*exp(-0.5*(xVals-mean_hist).^2/std_hist);
% %plot(xVals, f_norm, 'k--')
% xlabel('Innovation North [m/s]')
% legend('Simulation 1','Simulation 2')
% grid on;
% 
% subplot(2,1,2)
% ii = 1; dimen = 2;
% jj = 2;
% hist([squeeze(innovation(dimen,1:end-2,ii)); ...
%       squeeze(innovation(dimen,1:end-2,jj))]'); hold on;
% std_hist = std([squeeze(innovation(dimen,1:end-2,ii)), ...
%       squeeze(innovation(dimen,1:end-2,jj))]);
% mean_hist = mean([squeeze(innovation(dimen,1:end-2,ii)), ...
%       squeeze(innovation(dimen,1:end-2,jj))]);
% N_vals = 100;
% histMax = 50;
% xLims = xlim;
% xVals = linspace(xLims(1),xLims(2), N_vals);
% f_norm = histMax*1/(std_hist*sqrt(2*pi)).*exp(-0.5*(xVals-mean_hist).^2/std_hist);
% %plot(xVals, f_norm, 'k--')
% xlabel('Innovation East [m/s]')
% grid on;
% 
% print(strcat('fig/histogramm_innovation_cut_',figName),'-depsc')
% 
% % %% Rotational plot
% % N_t = length(time);
% % figure('Position',[100,200,800,700]);
% % set(groot,'DefaultAxesFontSize',14)
% % set(groot,'DefaultLineLineWidth',1.0)
% % subplot(2,1,1)
% % iter = 1; dim = 1;
% % plot(time, x_real(dim,1:N_t),'b-x'); hold on;
% % plot(x_simu(dim,1:N_t,iter),'ro'); hold on;
% % plot(time, x_KF(dim,1:N_t, iter),'r')
% % time_long = linspace(0,time(end),size(x_tild,2));
% % plot(time_long, x_tild (dim,1:end, iter),'r--')
% % subplot(2,1,2)
% % ylabel('North [m]'); 
% % iter = 1; dim = 2;
% % plot(time, x_real(dim,1:N_t),'b-x'); hold on;
% % plot(x_simu(dim,1:N_t,iter),'ro'); hold on;
% % plot(time, x_KF(dim,1:N_t, iter),'r')
% % time_long = linspace(0,time(end),size(x_tild,2));
% % plot(time_long, x_tild (dim,1:end, iter),'r--')
% % 
% % 
% % ylabel('East [m]'), xlabel('Time [s]')
% % legend('Real measurement', 'Simulation 1', 'Kalman Filter 1', 'Inner loop 1')
% % 
% % print(strcat('fig/kalmanFilter_sigma',figName),'-depsc')
% % %%
% % iter = 1; dim = 2;
% % 
% % figure;
% % plot(time,squeeze(innovation(dim,1:end,iter))); hold on;
% % plot([time(1), time(end)],[2*pi,2*pi],'k--')
% % plot([time(1), time(end)],(-1)*[2*pi,2*pi],'k--')
% end