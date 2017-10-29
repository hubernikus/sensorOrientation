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

%%
% Numerical values
omega = pi/100;     % [Hz]
r_circ = 500;       % [m]

vel0 = 0;
omega0 = omega;

dT = 0.2; % [s] - Time step

N_data = 1000; % number of generated data samples

%% Task 1 - Simulation of nominal measurment
close all; 

phi_dot = omega*dT*1:N_data;
x_inert = r_circ*[cos(phi_dot); sin(phi_dot)];

acc_body = [0;r_circ*omega0]*ones(1,N_data);
gyro_mb_body = omega0 *ones(1,N_data);


figure;
subplot(2,2,1);
plot(x(1,:),x(2,:),'color',[0.2,0.2,1])
xlabel('$e_1$ [m]','Interpreter','latex'); ylabel('$e_2$ [m]','Interpreter','latex')


%% Task 2 - 
