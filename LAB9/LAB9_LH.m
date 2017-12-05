%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Sensor Orientation - EPFL
%           LAB 9
%
%           Author: Huber Lukas
%           Date: 2017-11-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all; clear variables; clc;

addpath(genpath('functions'))

%% Rescaling sensor error
clc;
g = 9.81; % [m/s^2] - Gravity

%% Simulation Parameters
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
[acc, gyro, x_real, v_real, phi_real, time] = ... 
               IMUsens_simulation(samplingFreq, omega0, r_circ, phi0, azim_init);

% GPS Simulation 
simga_x_GPS = 0.5; simga_y_GPS = 0.5; % [m]
sigma_GPS = [simga_x_GPS,simga_y_GPS];
x_simu = x_real + [whiteNoiseGen(size(x_real,2),1,simga_x_GPS); whiteNoiseGen(size(x_real,2),1,simga_y_GPS)];

N_sim = size(x_simu,2);

%% Kalman Filtering
for iter = 1:2
    dt_kf = dt_kf_list(iter);
    
    % Kalman Filtering
    modelType = 'constAcc';
    [model, x_filt] = ...
            init_kalmanFilter(modelType, dt_kf, r_circ, omega0, sigma_GPS);
    
    [x_KF(:,:,iter), innovation(:,:,iter), sigma_pred(:,iter), x_tild] ...
                = kalmanFilter_extended(model, x_simu, x_filt, dt_kf, dt_gps);
end

sigma_at = model.sigma_at;

%%


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
axis equal;
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
for iter = 1:size(sigma_pred,2)
    plot(dt_gps*(0:size(sigma_pred,1)-1),sigma_pred(:,iter)); hold on;
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
hist([squeeze(innovation(dimen,1:end-2,ii)); ...
      squeeze(innovation(dimen,1:end-2,jj))]'); hold on;
std_hist = std([squeeze(innovation(dimen,1:end-2,ii)), ...
      squeeze(innovation(dimen,1:end-2,jj))]);
mean_hist = mean([squeeze(innovation(dimen,1:end-2,ii)), ...
      squeeze(innovation(dimen,1:end-2,jj))]);
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
hist([squeeze(innovation(dimen,1:end-2,ii)); ...
      squeeze(innovation(dimen,1:end-2,jj))]'); hold on;
std_hist = std([squeeze(innovation(dimen,1:end-2,ii)), ...
      squeeze(innovation(dimen,1:end-2,jj))]);
mean_hist = mean([squeeze(innovation(dimen,1:end-2,ii)), ...
      squeeze(innovation(dimen,1:end-2,jj))]);
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
hist([squeeze(innovation(dimen,1:end-2,ii)); ...
      squeeze(innovation(dimen,1:end-2,jj))]'); hold on;
std_hist = std([squeeze(innovation(dimen,1:end-2,ii)), ...
      squeeze(innovation(dimen,1:end-2,jj))]);
mean_hist = mean([squeeze(innovation(dimen,1:end-2,ii)), ...
      squeeze(innovation(dimen,1:end-2,jj))]);
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
hist([squeeze(innovation(dimen,1:end-2,ii)); ...
      squeeze(innovation(dimen,1:end-2,jj))]'); hold on;
std_hist = std([squeeze(innovation(dimen,1:end-2,ii)), ...
      squeeze(innovation(dimen,1:end-2,jj))]);
mean_hist = mean([squeeze(innovation(dimen,1:end-2,ii)), ...
      squeeze(innovation(dimen,1:end-2,jj))]);
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
% figure('Position',[100,200,800,600]);
% set(groot,'DefaultAxesFontSize',14)
% set(groot,'DefaultLineLineWidth',1.2)
% subplot(2,1,1)
% for iter = 1:size(innovation,2)
%     plot(dt_gps*(0:size(innovation,3)-1),squeeze(innovation(1,iter,:))); hold on;
% end
% ylabel('Innovation (east) [m]')
% legend('Prediction interval: 1 Hz','Prediction interval: 0.1 Hz')
% subplot(2,1,2)
% for iter = 1:size(innovation,2)
%     plot(dt_gps*(0:size(innovation,3)-1),squeeze(innovation(1,iter,:))); hold on;
% end
% ylabel('Innovation (north) [m]')
% xlabel('Time [s]')
% print('fig/allInovations','-depsc')

% %%
% figure('Position',[100,200,800,600]);
% set(groot,'DefaultAxesFontSize',14)
% set(groot,'DefaultLineLineWidth',1.2)
% subplot(2,1,1)
% for iter = 1:size(innovation,2)
%     sabilisationTime(iter) = find(abs((sigma_pred(iter,(2:end))-sigma_pred(iter,(1:end-1)))...
%                         ./sigma_pred(iter,(2:end))) < 1e-4 ,1);
%     plot(dt_gps*(sabilisationTime(iter)-1:size(innovation,3)-1),squeeze(innovation(1,iter,sabilisationTime(iter):end))); hold on;
% end
% ylabel('Innovation (east) [m]')
% legend('Prediction interval: 1 Hz','Prediction interval: 0.1 Hz')
% subplot(2,1,2)
% for iter = 1:size(innovation,2)
%     plot(dt_gps*(sabilisationTime(iter)-1:size(innovation,3)-1),squeeze(innovation(2,iter,sabilisationTime(iter):end))); hold on;
% end
% ylabel('Innovation (north) [m]')
% xlabel('Time [s]')
% print('fig/stableInovations','-depsc')

