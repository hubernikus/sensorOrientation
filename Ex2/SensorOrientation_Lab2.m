%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               Sensor Orientation - EPFL
%               Lab 2 - 
%               Lukas Huber - 2017/10/13
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accelerometer in x & z direction
% Sensors used: 3 - FSAS/IMAR (tactical, FOG) 6 - Smart-phone (low-cost,
% MEMS) 

clc; clear variables; close all;

addpath(genpath('../Ex1'))
addpath(genpath('../functions'))

dataPath = 'imuData2017/'; 

%% Define General Parameters 

c = ['r','g','b']; 


%% Notes Black-Board
%
% d = [t, x_0 y_0 z_0, x_A y_A z_A]
% x_k = exp(-beta Dt) * x_{k-1} + w_k (NOT \simga_w ?)
% (simga_GM)^2 = (exp(-beta Dt)^2 * \sigma_GM - (sigma_W)^2
% sigma_W^2 = sigma _{GM1}^2 - exp(-2 beta Dt) sigma_GM1^2
% simga_w = sqrt(1-exp(-2 beta Dt) sigma_GM^2
% w_k = sqrt(sigma_W^2) * randn 


%% Import Data
[dataIMAR, fIMU] = readimu(strcat(dataPath, 'imu3_20171006_imar-fsas.imu'),'IMAR');
load(strcat(dataPath, 'imu6_Phone.mat'));

%% IMAR Analysis
%N_sample = 'end';
N_max = 10^10;
N_sample = min(size(dataIMAR,1), N_max);

% Gyroscope Data
% figure(1001);
% for ii = 1:3
%     plot(dataIMAR(1:N_sample,1), dataIMAR(1:N_sample,ii+1),c(ii)); hold on;
% end
% xlim([dataIMAR(1,1),dataIMAR(N_sample,1)])

figure(1);
for ii = [1,3]
    plot(dataIMAR(1:N_sample,1)-dataIMAR(1,1), dataIMAR(1:N_sample,ii+4),c(ii)); hold on;
end
xlim([dataIMAR(1,1),dataIMAR(N_sample,1)]-dataIMAR(1,1))
grid on;
xlabel('Sample []','Interpreter','latex'); ylabel('Accelerometer measurement [$m/s^2$]','Interpreter','latex')
legend('x direction','z direction')

%% Phone Analysis
N_sample = min(size(t_a,1), N_max);

figure(2);
for ii = [1,3]
    plot(t_a(1:N_sample), a(1:N_sample,ii),c(ii)); hold on;
end
xlim([t_a(1,1),t_a(N_sample,1)])
xlabel('Sample []','Interpreter','latex'); 
ylabel('Accelerometer measurement [$m/s^2$]','Interpreter','latex')
legend('x direction','z direction')
grid on;


%% Correlation analysis


