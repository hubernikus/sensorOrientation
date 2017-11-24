%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Sensor Orientation - LAB 5 
%           EPFL
%
%           Author: Huber Lukas
%           Date: 2017-11-03
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all; clear all; clc;

addpath(genpath('../LAB1'))
addpath(genpath('../LAB3'))

%% Rescaling sensor errorc
clc;
g = 9.81; % [m/s^2] - Gravity

samplingFreq = 10;     % [Hz]
<<<<<<< HEAD
%samplingFreq = 100;    % [Hz]
=======
>>>>>>> 654334958b97321bbb0420ea5429a16f2dd048db

b_g = 10;               % [deg/h]
simga_G_GM = 0.005;     % [deg/s/sqrt(Hz)]
beta_G_inv = 100;       % [s]
sigma_G_vel = 0.1;      % [deg/sqrt(h)]
b_A = 1;                %[mg]
sigma_A_vel = 50;       % [mu g/ /sqrt(Hz)]

% Rescale to SI units
b_g = b_g/180*pi/3600 
simga_G_GM = simga_G_GM/180*pi*sqrt(samplingFreq)
beta_G_inv = beta_G_inv
sigma_G_vel = sigma_G_vel/180*pi*1/sqrt(3600)*sqrt(samplingFreq)
b_A = b_A * g * 1e-3
sigma_A_vel = sigma_A_vel * 1e-6 * g * sqrt(samplingFreq)

%% Ex 2 
% Initial conditions
omega0 = pi/100;     % [rad/s]
r_circ = 500;       % [m]
azim_init = 90/180*pi; % [rad]
phi0 = 90/180*pi;
x0 = [0;r_circ];
vel0 = [-omega0*r_circ;0]; % initial velocity
initHeading = [phi0+azim_init];

% Simulate Measurement
[acc, gyro, x_real, v_real, phi_real, time] = ... 
           IMUsens_simulation(samplingFreq, omega0, r_circ, phi0, azim_init);
             
% Add different noise on measurement
% Number of samples
N_samples = length(acc);

% Gyro bias
gyro_bias = ones(1, N_samples) * b_g;

% Gyro correlated noise
deltaT = 1/samplingFreq;
gyro_corrNoise =  gaussMarkovGen(N_samples, 1, deltaT, beta_G_inv, simga_G_GM);

% Gyro random walk
gyro_randomWalk = whiteNoiseGen(N_samples, 1, sigma_G_vel);

gyro_errors = gyro + gyro_bias + gyro_corrNoise + gyro_randomWalk;

% Accelerometer bias
acc_bias = ones(2, N_samples) * b_A;

% Accelerometer noise
acc_noise = whiteNoiseGen(N_samples, 2, sigma_A_vel);

acc_errors = acc + acc_bias + acc_noise; 

%%
close all;

i = 1;
intOrder = 2;
[x_sim{i}, v_sim{i}, phi_sim{i}] = inertialNavigation(x0, vel0, initHeading, acc, gyro_errors, time,intOrder);

i = 2; 
[x_sim{i}, v_sim{i}, phi_sim{i}] = inertialNavigation(x0, vel0, initHeading, acc_errors, gyro, time,intOrder);

i = 3
[x_sim{i}, v_sim{i}, phi_sim{i}] = inertialNavigation(x0, vel0, initHeading, acc, gyro, time,intOrder);


% Calculate Error - Postition
posErr = []; velErr = []; angErr=[];

i = 3;
titleName = 'no_errors';
[posErr(i), velErr(i), angErr(i)] = createErrorPlot(x_sim{i}, x_real, v_sim{i}, v_real, phi_sim{i}, phi_real, titleName)

i = 1;
titleName = 'Gyro_Errors';
[posErr(i), velErr(i), angErr(i)] = createErrorPlot(x_sim{i}, x_real, v_sim{i}, v_real, phi_sim{i}, phi_real, titleName)

i = 2;
titleName = 'acc_Errors';
[posErr(i), velErr(i), angErr(i)] = createErrorPlot(x_sim{i}, x_real, v_sim{i}, v_real, phi_sim{i}, phi_real, titleName)


%%
fprintf('gyro_bias: \n')
i = 4;
gyro_partError = gyro + gyro_bias;
[x_sim, v_sim, phi_sim] = inertialNavigation(x0, vel0, initHeading, acc, gyro_partError, time,intOrder);
[posErr(i), velErr(i), angErr(i)] = createErrorPlot(x_sim, x_real, v_sim, v_real, phi_sim, phi_real)

fprintf('gyro_corrNoise: \n')
i = 5;
gyro_partError = gyro + gyro_corrNoise;
[x_sim, v_sim, phi_sim] = inertialNavigation(x0, vel0, initHeading, acc, gyro_partError, time,intOrder);
[posErr(i), velErr(i), angErr(i)] = createErrorPlot(x_sim, x_real, v_sim, v_real, phi_sim, phi_real)

fprintf('gyro_randomWalk: \n')
i = 6;
gyro_partError = gyro + gyro_randomWalk;
[x_sim, v_sim, phi_sim] = inertialNavigation(x0, vel0, initHeading, acc, gyro_partError, time,intOrder);
[posErr(i), velErr(i), angErr(i)] = createErrorPlot(x_sim, x_real, v_sim, v_real, phi_sim, phi_real)

fprintf('acc_bias: \n')
i = 7;
acc_partError = acc + acc_bias; 
[x_sim, v_sim, phi_sim] = inertialNavigation(x0, vel0, initHeading, acc_partError, gyro, time,intOrder);
[posErr(i), velErr(i), angErr(i)] = createErrorPlot(x_sim, x_real, v_sim, v_real, phi_sim, phi_real)

fprintf('acc_noise: \n')
i = 8;
acc_partError = acc + acc_noise; 
[x_sim, v_sim, phi_sim] = inertialNavigation(x0, vel0, initHeading, acc_partError, gyro, time,intOrder);
[posErr(i), velErr(i), angErr(i)] = createErrorPlot(x_sim, x_real, v_sim, v_real, phi_sim, phi_real)


%% Print table to LATEX

fileID = fopen('table_errors.tex','w');
fprintf(fileID,'Position Errors [m]& %3.2f & %3.2f & %3.2f & %3.2f & %3.2f & %3.2f  & %3.2f  & %3.2f  \\\\ \\hline \n',posErr);
fprintf(fileID,'Velocity Errors [m]& %3.2f & %3.2f & %3.2f & %3.2f & %3.2f & %3.2f  & %3.2f  & %3.2f  \\\\ \\hline \n',velErr);
fprintf(fileID,'Azimuth Error [m]& %3.2f & %3.2f & %3.2f & %3.2f & %3.2f & %3.2f  & %3.2f  & %3.2f  \\\\ \\hline ',angErr);
fclose(fileID);