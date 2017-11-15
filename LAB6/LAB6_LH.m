%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Sensor Orientation - EPFL
%           LAB 6 
%
%           Author: Huber Lukas
%           Date: 2017-11-10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all; clear variables; clc;

addpath(genpath('../functions'))
%addpath(genpath('../LAB1'))
%addpath(genpath('../LAB3'))

% Import Data
[dataIXSEA, fIXSEA] = readimu('data/ixsea_2nd_group.imu','IXSEA');

% Initial data
g = 9.81; % [m/s^2] - Gravity

% Start and stop time in TOF (Time of the week)
timeLimit = [480610,480710];

% Initial Orientation
heading  = 299.337; % [deg]
roll = 2.524; % [deg]
pitch = 1.710; % [deg]

%%
% After observing the data, measurements was limited between [10 - 80s]
t_smooth = [10,80];

% Find limit times
indLimit = [find(dataIXSEA(:,1) > timeLimit(1) + t_smooth(1), 1), ...
            find(dataIXSEA(:,1) > timeLimit(1) + t_smooth(2), 1)];

% Time set to zero at timeLimit(1) = 480610
time = dataIXSEA(indLimit(1):indLimit(2),1)-(timeLimit(1)++ t_smooth(1)); 

gyro = dataIXSEA(indLimit(1):indLimit(2),2:4);
acc = dataIXSEA(indLimit(1):indLimit(2),5:7);

% Demean data 
accMean = mean(acc);
gyroMean = mean(gyro);

acc = acc - accMean;
gyro = gyro - gyroMean;

% Place Data in NED Frame
acc(:,2) = -acc(:,2); acc(:,3) = -acc(:,3); % NWU - NED
gyro(:,2) = -gyro(:,2); % NWD - NED

% Sampling Frequency
dTime = (time(2:end)-time(1:end-1));     % [Hz]
stdTime = std(dTime);
T_sampling = mean(dTime);

if(stdTime > T_sampling *1e-3)
    warning('The data has highly varying sampling frequency. \n')
end

%% Plot Signal

close all;

% Plot results
if(0) % 1 to run and save again
    figure('Position',[0,0,800,800]); % Plot measurements
    set(groot,'DefaultAxesFontSize',14)
    set(groot,'DefaultLineLineWidth',1.2)
    
    h1 = subplot(2,1,1);
    plot(h1, time,acc(:,1)-mean(acc(:,1))); hold on;
    plot(h1, time,acc(:,2)-mean(acc(:,2))); hold on;
    plot(h1, time,acc(:,3)-mean(acc(:,3))); hold on;
    xlim([time(1), time(end)]);
    xlabel('Time [s]'), ylabel('Accelerometer measurement [m/s^2]')
    legend('x direction','y direction','z direction')
    
    h2 = subplot(2,1,2);
    plot(h2, time,gyro(:,1)-mean(gyro(:,1))); hold on;
    plot(h2, time,gyro(:,2)-mean(gyro(:,2))); hold on;
    plot(h2, time,gyro(:,3)-mean(gyro(:,3))); hold on;
    xlim([time(1), time(end)]);
    xlabel('Time [s]'), ylabel('Gyroscope measurement [rad/s]')
    
    print('fig/IMU_measurement','-depsc')
end

%% Calculate norm of the signals
accNorm = norm(accMean)
gyroNorm = norm(gyroMean)

gravityTheoretical = -(980000 + 550)* 1e-5; %[m/s^2] 
angularRotation = 7.2921150 * 1e-5; %[rad/s] 

varAcc = (accNorm - abs(gravityTheoretical))/abs(gravityTheoretical)
varGyr = (gyroNorm -angularRotation )/angularRotation




%% Important

% Take mean before norm, for average signal

%% Important




%% Ex 2 
% Initial conditions
omega0 = pi/100;     % [rad/s]
r_circ = 500;       % [m]
azim_init = 90/180*pi; % [rad]
phi0 = 90/180*pi;
x0 = [0;0;0];
vel0 = [-omega0*r_circ;0]; % initial velocity
initHeading = [phi0+azim_init];

            
% Add different noise on measurement
% Number of samples
%N_samples = length(acc);



%%
close all;
i = 1; 
[x_sim{i}, v_sim{i}, phi_sim{i}] = inertialNavigation(x0, vel0, initHeading, acc, gyro, time, intOrder);

titleName = 'test';
[posErr(i), velErr(i), angErr(i)] = createErrorPlot(x_sim{i}, x_real, v_sim{i}, v_real, phi_sim{i}, phi_real, titleName)



%% Print table to LATEX

% fileID = fopen('table_errors.tex','w');
% fprintf(fileID,'Position Errors [m]& %3.2f & %3.2f & %3.2f & %3.2f & %3.2f & %3.2f  & %3.2f  & %3.2f  \\\\ \\hline \n',posErr);
% fprintf(fileID,'Velocity Errors [m]& %3.2f & %3.2f & %3.2f & %3.2f & %3.2f & %3.2f  & %3.2f  & %3.2f  \\\\ \\hline \n',velErr);
% fprintf(fileID,'Azimuth Error [m]& %3.2f & %3.2f & %3.2f & %3.2f & %3.2f & %3.2f  & %3.2f  & %3.2f  \\\\ \\hline ',angErr);
% fclose(fileID);