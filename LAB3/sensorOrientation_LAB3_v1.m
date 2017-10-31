%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Sensor Orientation -- EPFL 
%                   LAB 4 
%                       
%                   Author: Lukas Huber
%                   Date: 2017/10/27
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear variables;

%% Simulation

% Initial conditions
omega0 = pi/100;     % [rad/s]
r_circ = 500;       % [m]
azim_init = 90/180*pi; % [rad]
phi0 = 90/180*pi;
x0 = [0;r_circ];
vel0 = [-omega0*r_circ;0]; % initial velocity
initHeading = [phi0+azim_init];

% Create variables
x_sim =[]; v_sim =[]; phi_sim = [];

% Create measurement 1 - 10 Hz
sampleRate = 10; % [Hz]
[acc_body_10hz, gyro_mb_10hz, x_real_10hz, v_real_10hz, phi_real_10hz, time_10hz] = IMUsens_simulation(sampleRate, omega0, r_circ, phi0, azim_init);

intOrder = 1;
i = 1;
[x_sim{i}, v_sim{i}, phi_sim{i}] = inertialNavigation(x0, vel0, initHeading, acc_body_10hz, gyro_mb_10hz, time_10hz,intOrder);

intOrder = 2;
i = 2;
[x_sim{i}, v_sim{i}, phi_sim{i}] = inertialNavigation(x0, vel0, initHeading, acc_body_10hz, gyro_mb_10hz, time_10hz,intOrder);

% Create measurement 2 - 100Hz
sampleRate100 = 100;
[acc_body_100hz, gyro_mb_100hz, x_real_100hz, v_real_100hz, phi_real_100hz, time_100hz] = IMUsens_simulation(sampleRate100, omega0, r_circ, phi0, azim_init);

intOrder = 1;
i=3;
[x_sim{i}, v_sim{i}, phi_sim{i}] = inertialNavigation(x0, vel0, initHeading, acc_body_100hz, gyro_mb_100hz, time_100hz,intOrder);

intOrder = 2;
i=4;
[x_sim{i}, v_sim{i}, phi_sim{i}] = inertialNavigation(x0, vel0, initHeading, acc_body_100hz, gyro_mb_100hz, time_100hz,intOrder);


%% Calculate Error - Postition

lineName = {'10Hz, 1st order', '10Hz, 2st order','100Hz, 1st order','100Hz, 2st order'};
close all;

max_posErrs = zeros(length(x_sim),3);

figure('Position',[0,0,800,1200]); % Plot results
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)

err_pos = [];
h1 = subplot(3,1,1);
ylabel('Position error east [m]')
h2 = subplot(3,1,2);
ylabel('Position error north [m]')
h3 = subplot(3,1,3);
ylabel('Total position Error [m]'); xlabel('Azimuth [deg]')

for i = 1:length(x_sim)
    if i <=2 % 10hz
        x_real = x_real_10hz;
        phi_real = phi_real_10hz*180/pi;
    else
        x_real = x_real_100hz;
        phi_real = phi_real_100hz*180/pi;
    end
          
    err_pos{i} = x_sim{i}-x_real;
    hold(h1,'on')
    plot(h1, phi_real, err_pos{i}(1,:)); 
    hold(h2,'on')
    plot(h2, phi_real, err_pos{i}(2,:));  
    hold(h3,'on')
    absErr = sqrt(sum(err_pos{i}.^2,1));
    plot(h3, phi_real, absErr); 
    
    % Maximal errors
    max_posErrors(i,1) = max(abs(err_pos{i}(1,:)));
    max_posErrors(i,2) = max(abs(err_pos{i}(2,:)));
    max_posErrors(i,3) = max(absErr);
end

hold off;
subplot(3,1,2)
legend(lineName,'Location','northwest')

print(strcat('fig/','Error_position'),'-depsc')


%% Calculate Error -- Velocity

figure('Position',[0,0,800,1200]); % Plot results
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)

max_velErrors = zeros(length(x_sim),3);

error_vel = [];
h1 = subplot(3,1,1);
ylabel('Velocity error east [m/s]')
h2 = subplot(3,1,2);
ylabel('Veloctiy erro north [m/s]')
h3 = subplot(3,1,3);
ylabel('Total velocty Error [m/s]'); xlabel('Azimuth [deg]')

for i = 1:length(x_sim)
    if i <=2 % 10hz
        v_real = v_real_10hz;
        phi_real = phi_real_10hz*180/pi;
    else
        v_real = v_real_100hz;
        phi_real = phi_real_100hz*180/pi;
    end
          
    err_vel{i} = v_sim{i}(:,1:end-1)-v_real;
    hold(h1,'on')
    plot(h1, phi_real(1:end-1), err_vel{i}(1,:)); 
    hold(h2,'on')
    plot(h2, phi_real(1:end-1), err_vel{i}(2,:));  
    hold(h3,'on')
    absErr = sqrt(sum(err_vel{i}.^2,1));
    plot(h3, phi_real(1:end-1), absErr); 
    
    % Maximal errors
    max_velErrors(i,1) = max(abs(err_vel{i}(1,:)));
    max_velErrors(i,2) = max(abs(err_vel{i}(2,:)));
    max_velErrors(i,3) = max(absErr);
end

hold off;
subplot(3,1,2)
legend(lineName,'Location','west')

print(strcat('fig/','Error_vel'),'-depsc')


%% Calculate Error -- Azimuth

figure('Position',[0,0,800,600]);  % Plot results
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)

max_angErrors = zeros(length(x_sim),1);

error_vel = [];

for i = 1:length(x_sim)
    if i <=2 % 10hz
        phi_real = phi_real_10hz*180/pi;
    else
        phi_real = phi_real_100hz*180/pi;
    end
    
    err_phi{i} = phi_sim{i}*180/pi-phi_real;
    absErr = abs(err_phi{i}); 
    max_angErrors(i) = max(absErr);
    
    plot(phi_real, absErr); hold on; 
    
end

ylabel('Total azimuth error [deg]'); xlabel('Real Azimuth [deg]')
legend(lineName,'Location','west')

print(strcat('fig/','Error_ang'),'-depsc')


%% Print table to LATEX

fileID = fopen('table_errors.tex','w');
fprintf(fileID,'Error x [m]& %3.2f & %3.2f & %3.2f & %3.2f  \\\\ \\hline \n',max_posErrors(:,1));
fprintf(fileID,'Error y [m]& %3.2f & %3.2f & %3.2f & %3.2f  \\\\ \\hline \n',max_posErrors(:,2));
fprintf(fileID,'Error total [m]& %3.2f & %3.2f & %3.2f & %3.2f  \\\\ \\hline  \\hline  \n',max_posErrors(:,3));
pow = -5;
fprintf(fileID,'Error $v_x $ [1e%d m/s] & %3.2f & %3.2f& %3.2f & %3.2f\\\\ \\hline \n',pow, max_velErrors(:,1)*10^(-pow));
pow = -2;
fprintf(fileID,'Error  $v_y$ [1e%d m/s] & %3.2f & %3.2f& %3.2f & %3.2f\\\\ \\hline \n',pow, max_velErrors(:,2)*10^(-pow));
fprintf(fileID,'Error $v_{tot} $ [1e%d m/s] & %3.2f & %3.2f& %3.2f & %3.2f\\\\ \\hline \\hline \n',pow, max_velErrors(:,3)*10^(-pow));
pow = -10;
fprintf(fileID,'Angle [1e%d deg]& %3.2f & %3.2f & %3.2f & %3.2f \\\\ \\hline',pow, max_angErrors*10^(-pow));
fclose(fileID);