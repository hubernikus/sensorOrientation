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

addpath(genpath('../functions'))

%%
% Colors used for plots
cols = [1 0 0; 0 1 0; 0 0 1; 0.7 1 0.7; 0.5 0.7 1];

% Initial conditions
phi0 = 0; 
omega0 = pi/100;     % [Hz]
r_circ = 500;       % [m]
azim_init = 90*pi/180; % [rad]
x0 = [r_circ; 0];
vel0 = [0; omega0*r_circ]; % initial velocity


% Simulation paramters
dT = 0.2; % [s] - Time step
N_data = 2000; % number of generated data samples


% Task 1 - Simulation of nominal measurment
[acc_body, gyro_mb] = IMUsens_simulation(N_data, dT, oemga0, r_circ);


% Task 2 - tratdown inertial navigation
[] = inertialNavigation(x0, vel0, acc_body, gyro_mb, dT)

omega = omega0* ones(1,N_data);
t = dT*(1:N_data+1);

delta_phi = omega.*(t(2:end)-t(1:end-1))+phi0;
phi = cumsum(delta_phi);
x_m = r_circ*[cos(phi_dot); sin(phi_dot)];

%% 
close all;

figure;
subplot(2,2,1);
plot(x_m(1,:),x_m(2,:),'color',[0.2,0.2,1])
xlabel('$e_1$ [m]','Interpreter','latex'); ylabel('$e_2$ [m]','Interpreter','latex')


figure;
plot(t,x_m(1,:),'color',cols(1,:)); hold on;
plot(t,x_m(2,:),'color',cols(2,:))
plot(t,acc_body(1,:),'color',cols(3,:))
plot(t,acc_body(2,:),'color',cols(4,:))
plot(t,gyro_mb,'color',cols(5,:))
xlabel('Time [s]')
legend({'Pos $x_1$','Pos $x_1$','Accelerometer $x_1$','Accelerometer $x_2$','Gyroscope'},'Interpreter','latex')

%% Task 2 - 
