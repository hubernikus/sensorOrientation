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

% Import Data
[dataIXSEA, fIXSEA] = readimu('data/ixsea_2nd_group.imu','IXSEA');

% Start and stop time in TOF (Time of the week)
timeLimit = [480610,480710];

% Rotation Matrix of Sensor Reading
R_NWU_NED = [1 0 0; 0 -1 0; 0 0 -1];

% Initial Orientation
roll_ref = 2.524*pi/180; % [rad] 
pitch_ref = 1.710*pi/180*(-1); % [rad]-- NWD -> NED frame
heading_ref  = 299.337*pi/180; % [rad]

% Rotation Matrix of Reference Inputs
R_NWD_NED = [1 0 0; 0 -1 0; 0 0 1];

% Lattitue EPFL
% N 46°31'17''
phi_ref = pi/180*(46 + 1/60*(31 + 1/60*17) );


%%
% After observing the data, measurements was limited between [10 - 80s]
t_smooth = [10,80];

% Find limit times
indLimit = [find(dataIXSEA(:,1) > timeLimit(1) + t_smooth(1), 1), ...
            find(dataIXSEA(:,1) > timeLimit(1) + t_smooth(2), 1)];

% Time set to zero at timeLimit(1) = 480610
time = dataIXSEA(indLimit(1):indLimit(2),1)-(timeLimit(1)++ t_smooth(1)); 

gyro = R_NWU_NED*dataIXSEA(indLimit(1):indLimit(2),2:4)';
acc = R_NWU_NED*dataIXSEA(indLimit(1):indLimit(2),5:7)';

% Demean data 
accMean = mean(acc,2);
gyroMean = mean(gyro,2);

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
    plot(h1, time,acc(1,:)-mean(acc(1,:))); hold on;
    plot(h1, time,acc(2,:)-mean(acc(2,:))); hold on;
    plot(h1, time,acc(3,:)-mean(acc(3,:))); hold on;
    xlim([time(1), time(end)]);
    xlabel('Time [s]'), ylabel('Accelerometer measurement [m/s^2]')
    legend('x direction','y direction','z direction')
    
    h2 = subplot(2,1,2);
    plot(h2, time,gyro(1,:)-mean(gyro(1,:))); hold on;
    plot(h2, time,gyro(2,:)-mean(gyro(2,:))); hold on;
    plot(h2, time,gyro(3,:)-mean(gyro(3,:))); hold on;
    xlim([time(1), time(end)]);
    xlabel('Time [s]'), ylabel('Gyroscope measurement [rad/s]')
    
    print('fig/IMU_measurement','-depsc')
end

%% Calculate norm of the signals
% Take mean before norm, for average signal

accNorm = norm(accMean);
gyroNorm = norm(gyroMean);

gravityTheoretical = -(980000 + 550)* 1e-5; %[m/s^2] 
angularRotation = 7.2921150 * 1e-5; %[rad/s] 

% Percental differenc
varAcc = (accNorm - abs(gravityTheoretical))/abs(gravityTheoretical);
varGyr = (gyroNorm -angularRotation )/angularRotation;

fprintf('\n')
fprintf('Accelerometer Measuerement: %0.7e deg, Relavitve error: %0.4e \n',accNorm, varAcc)
fprintf('Accelerometer Bias: %0.7e m/s^2 \n', abs(accNorm-abs(gravityTheoretical)))
fprintf('Gyroscope Masurement: %0.7e deg, Relavitve error: %0.4e \n',gyroNorm, varGyr)
fprintf('Gyroscope Bias: %0.7e rad/s \n', abs(gyroNorm-angularRotation) )
fprintf('\n')

%% IV - Leveling Accelerometer to NED
% Important

% Pitch 
roll = asin(accMean(2)/-accNorm);

R_p = [1 0 0;
        0 cos(roll) sin(roll);
        0 -sin(roll) cos(roll)];

% Roll
pitch = asin(accMean(1)/accNorm);

R_r = [cos(pitch) 0 -sin(pitch);
        0 1 0;
       sin(pitch) 0 cos(pitch)];

% Compare relative Error
roll_relError = (roll-roll_ref)/roll_ref;
pitch_relError = (pitch-pitch_ref)/pitch_ref;

   
% R = R_r * R_p * R_y
accLeveled = (R_r * R_p)'*accMean;


fprintf('Roll: %3.4f deg, Reference: %3.4f, Relavitve error: %0.4e \n', ...
                        roll*180/pi, roll_ref*180/pi, roll_relError)
fprintf('Pitch: %3.4f deg, Reference: %3.4f, Relavitve error: %0.4e \n',...
                    pitch*180/pi, pitch_ref*180/pi, pitch_relError)
fprintf('\n')


%% IV - Gyrocompassing to Esimtate YAW
gyroLeveled = (R_r * R_p)'*gyroMean;


Az = atan2(-gyroLeveled(2),gyroLeveled(1))+2*pi;
R_y = [cos(Az) sin(Az) 0;
       -sin(Az) cos(Az) 0;
        0 0 1];

az_relError = (Az-heading_ref)/heading_ref;

fprintf('Azimuth: %3.4f deg, Reference: %3.4f, Relavitve error: %0.7e \n', ...
                        Az*180/pi, heading_ref*180/pi, az_relError)
fprintf('\n')

gyroCompassed = R_y' * gyroLeveled;

%% 10. Determine Lattitude Value
phi_1 = acos(gyroCompassed(1)/angularRotation);
phi_2 = asin(-gyroCompassed(3)/angularRotation);

fprintf('Ref: %2.7f, Lattitude 1: %2.7f, Lattitude 2: %2.7f \n' , ...
    phi_ref*180/pi, phi_1*180/pi, phi_2*180/pi)

fprintf('Difference -  Lattitude 1: %0.7e, Lattitude 2: %0.7e \n' , ...
     (phi_1-phi_ref), (phi_2-phi_ref))
