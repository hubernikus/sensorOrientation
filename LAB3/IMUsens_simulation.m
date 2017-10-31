function [acc, gyro, x_real, y_real] = IMUsens_simulation(sampleRate, omega, r_circ, phi0)
% Simulates the inertial sensors

% Sample time
%sampleRate = 1; % [Hz]

% Real position
t_final = 2*pi/omega;
time = 0:1/sampleRate:t_final;

x_real = r_circ*[cos(time*omega+phi0); sin(time*omega+phi0)];
y_real = [x_real_mb(:,2:end)-x_real_mb(:,1:end-1)]*sampleRate;
, r_circ)
% Simulates the inertial sensors
%dT = 1/sampleRate; % Time Step

N = length(time); % number of generated data samples

%dT = 1/sampleRate;

% for 2D
if(length(omega) == 1)
    omega = omega* ones(1,N);
end
acc = [zeros(1,length(omega));r_circ*omega.^2];
gyro = omega;

end