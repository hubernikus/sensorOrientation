function [acc, gyro] = IMUsens_simulation(N, sampleRate, omega, r_circ)
% Simulates the inertial sensors

dT = 1/sampleRate;

% for 2D
if(length(omega) == 1)
    omega = omega* ones(1,N);
end
acc = [zeros(1,length(omega));r_circ*omega.^2];
gyro = omega;

end