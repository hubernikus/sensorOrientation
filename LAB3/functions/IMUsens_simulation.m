function [acc, gyro] = IMUsens_simulation(N, dT, omega, r_circ)
% Simulates the inertial sensors

% for 2D
if(length(omega) == 1)
    omega = omega0* ones(1,N_data);
end
acc = [zeros(1,length(omega));r_circ*omega];
gyro = omega;

end