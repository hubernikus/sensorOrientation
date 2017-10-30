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

%addpath(genpath('functions'))

%%
% Colors used for plots
cols = [1 0 0; 0 1 0; 0 0 1; 0.7 1 0.7; 0.5 0.7 1];

% Initial conditions
phi0 = 0; 
omega0 = pi/100;     % [rad/s]
r_circ = 500;       % [m]
azim_init = 90/180*pi; % [rad]
phi0 = 0;
x0 = [0;r_circ];
vel0 = [-omega0*r_circ;0]; % initial velocity
initHeading = [phi0+azim_init];

% Sample time
sampleRate = 100; % [Hz]

% Real position
t_final = 2*pi/omega0;
time = 0:1/sampleRate:t_final;
x_real_mb = r_circ*[cos(time*omega0+initHeading); sin(time*omega0+initHeading)];
v_real_mb = [x_real_mb(:,2:end)-x_real_mb(:,1:end-1)]*sampleRate;

dT = 1/sampleRate; % Time Step

N_data = length(time); % number of generated data samples

% Task 1 - Simulation of nominal measurment
[acc_body_10hz, gyro_mb_10hz] = IMUsens_simulation(N_data, sampleRate, omega0, r_circ);

% Task 2 - tratdown inertial navigation
intOrder = 1;
[x_sim, v_sim, phi_sim] = inertialNavigation(x0, vel0, initHeading, acc_body_10hz, gyro_mb_10hz, time,intOrder);

%%
close all;
figure;
subplot(1,2,1)
plot(x_real_mb(1,:),x_real_mb(2,:),'b'); hold on;
plot(x_real_mb(1,1),x_real_mb(2,1),'bx')
plot(x_real_mb(1,end),x_real_mb(2,end),'bo')

plot(x_sim(1,:),x_sim(2,:),'r-'); 
plot(x_sim(1,1),x_sim(2,1),'rx')
plot(x_sim(1,end),x_sim(2,end),'ro')

xlabel('Position x [m]'); ylabel('Position z [m]');


subplot(1,2,2)
plot(v_real_mb(1,:),v_real_mb(2,:),'b'); hold on;
plot(v_real_mb(1,1),v_real_mb(2,1),'bx')
plot(v_real_mb(1,end),v_real_mb(2,end),'bo')

plot(v_sim(1,:),v_sim(2,:),'rx'); 
plot(v_sim(1,1),v_sim(2,1),'r-')
plot(v_sim(1,end),v_sim(2,end),'ro')

xlabel('Speed x [m/s]'); ylabel('Speed z [m/s]');





figure;
subplot(3,1,1)
plot(time,x_real_mb(1,:),'b'); hold on;
plot(time,x_real_mb(2,:),'r')
plot(time,x_sim(1,:),'b--');
plot(time,x_sim(2,:),'r--'); 
legend('Real x', 'Real y', 'Simulated x', 'Simulated y')
ylabel('Position')

subplot(3,1,2)
plot(time(1:end-1),v_real_mb(1,:),'b'); hold on;
plot(time(1:end-1),v_real_mb(2,:),'r')
plot(time(1:end),v_sim(1,:),'b--');
plot(time(1:end),v_sim(2,:),'r--'); 
legend('Real x', 'Real y', 'Simulated x', 'Simulated y')
ylabel('Velocity [m/s]')

subplot(3,1,3)
plot(time, time*omega0+initHeading,'b','linewidth',2); hold on;
plot(time,phi_sim(:),'r'); 
xlabel('Time[s]'); ylabel('Rotation [rad]');


%%
figure;

v2 = (x_sim(:,2:end)-x_sim(:,1:end-1))/dT;

plot(v2(1,:));hold on;
plot(v2(2,:)); 

%%
% 
% figure;
% 
% v2 = (x_real_mb(:,2:end)-x_real_mb(:,1:end-1))/dT;
% 
% plot(v2(1,:),v2(2,:),'rx'); hold on;
% 
% plot(v2(1,1),v2(2,1),'r-')
% plot(v2(1,end),v2(2,end),'ro')
% 
% xlabel('Speed x [m/s]'); ylabel('Speed z [m/s]');
% 
% % figure;
% plot(v_sim(1,:),v_sim(2,:))
% 
% figure;
% plot(time, phi_sim(1,:))

% %% 
% close all;
% 
% figure;
% subplot(2,2,1);
% plot(x_m(1,:),x_m(2,:),'color',[0.2,0.2,1])
% xlabel('$e_1$ [m]','Interpreter','latex'); ylabel('$e_2$ [m]','Interpreter','latex')
% 
% 
% figure;
% plot(t,x_m(1,:),'color',cols(1,:)); hold on;
% plot(t,x_m(2,:),'color',cols(2,:))
% plot(t,acc_body(1,:),'color',cols(3,:))
% plot(t,acc_body(2,:),'color',cols(4,:))
% plot(t,gyro_mb,'color',cols(5,:))
% xlabel('Time [s]')
% legend({'Pos $x_1$','Pos $x_1$','Accelerometer $x_1$','Accelerometer $x_2$','Gyroscope'},'Interpreter','latex')
% 
%% Task 2 - 
