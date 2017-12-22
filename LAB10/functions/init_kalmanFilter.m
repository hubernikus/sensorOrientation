function [model, X_init, dX_init] ...
            = init_kalmanFilter(dt_kf, r_circ, omega0, imu0, sigma_GPS, p_GPS_0, err)
% A priori Info Conditions
g = 9.81;
samplingFrequency = 1/dt_kf;

model = [];       
model.dimZ = 2;
model.dimX = 9;

% x [n x 1]- State vector
% x = [alpha, v_    n, v_e, p_n, p_e, ...]';
dX_init = [3*pi/180;-2;-1;0;0; ...
           err.gyro_bias; err.gyro_GM1.initSigma; ...
           err.acc1_GM1.initSigma; err.acc2_GM1.initSigma]; % initial error 
X_init = [90*pi/180; 0; omega0*r_circ; p_GPS_0(1); p_GPS_0(2); 0; 0; 0; 0];
X_init = X_init+dX_init;

alpha0 = X_init(1);

R_3 = [cos(alpha0), -sin(alpha0); sin(alpha0),cos(alpha0)];
acc_m = R_3*imu0(2:3)

F11 = [0 0 0 0 0; 
       -acc_m(1) 0 0 0 0;
       acc_m(2) 0 0 0 0;
       0 1 0 0 0;
       0 0 1 0 0];

F21 = zeros(4,5);

F12 = [1 1 0 0; 
       0 0 cos(alpha0) -sin(alpha0);
       0 0 sin(alpha0) cos(alpha0);
       0 0 0 0;
       0 0 0 0];

F22 = diag([0, -err.gyro_GM1.beta, -err.acc1_GM1.beta, -err.acc2_GM1.beta]);

model.F = [F11,F12; F21,F22];

G11 = [1 0 0;
       0 cos(alpha0) -sin(alpha0);
       0 sin(alpha0) cos(alpha0);
       0 0 0;
       0 0 0];
   
G22 = [0 0 0;
       1 0 0;
       0 1 0;
       0 0 1];

model.G = blkdiag(G11,G22);

% Uncertainty of motion model
model.W = diag([err.gyro_rw^2, err.acc1_GM1.sigma^2, err.acc2_GM1.sigma^2, ...
                2*err.gyro_GM1.sigma^2*err.gyro_GM1.beta, ...
                2*err.acc1_GM1.sigma^2*err.acc1_GM1.beta, ...
                2*err.acc2_GM1.sigma^2*err.acc2_GM1.beta]);


model.dimPhi = 9;
A = [-model.F, model.G*model.W*model.G'; zeros(model.dimPhi), model.F']*dt_kf;

B = expm(A);

model.Phi0 = B(model.dimPhi+1:end,model.dimPhi+1:end)';
model.Q0 = model.Phi0*B(model.dimPhi+1:end,1:model.dimPhi)';

% R [l x l] - Covariance matrix of measurment
model.R = diag([sigma_GPS(1)^2, sigma_GPS(2)^2]); % Diagonal elements negelcected

% z [l x 1] - Measurements 

% H [l x n]- Design Matrix ( z = H*x + v )
model.H = [0 0 0 1 0 0 0 0 0;
           0 0 0 0 1 0 0 0 0];

% P [n x n]- Covariance Matrix of the state vector
model.sigma.alpha = 2*pi/180;       % [rad]
model.sigma.v = 50;                 % [m/s]
model.sigma.p = 10;                 % [m]
model.sigma.b_c = 0.05*180/pi; % [rad/s]
model.sigma.b_g = 0.01*180/pi*sqrt(samplingFrequency);
model.sigma.b_A = 300 * 1e-6*g; % [m/s^2]

model.P0 = diag([model.sigma.alpha^2, ...
                 model.sigma.v^2, model.sigma.v^2,...
                 model.sigma.p^2, model.sigma.p^2, ...
                 model.sigma.b_c^2, model.sigma.b_g^2, ...
                 model.sigma.b_A^2, model.sigma.b_A^2 ]);
        
end