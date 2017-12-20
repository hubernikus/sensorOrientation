function [acc, gyro, x_real, v_real, phi_real, time] = ...
            IMUsens_simulation(sampleRate, omega0, r_circ, phi0, azim_init, t_final)
% Simulates the inertial sensors

% Real position
if nargin<6
    t_final = 2*pi/omega0;
end

time = 0:1/sampleRate:t_final+2;

% Simulation of x and v -- x_n, x_e
x_real = r_circ*[cos(time*omega0); sin(time*omega0)];
v_real = [x_real(:,2:end)-x_real(:,1:end-1)]*sampleRate;
x_real = x_real(:,1:end-1);

N = length(time); % number of generated data samples

% for 2D
if(length(omega0) == 1)
    omega = omega0* ones(1,N);
end
%acc = [-r_circ*omega.^2;zeros(1,length(omega))];
acc = [zeros(1,length(omega)); r_circ*omega.^2];
gyro = omega;

% Simulate Angle
initHeading = phi0+azim_init;
phi_real = time*omega0+initHeading-pi;

end