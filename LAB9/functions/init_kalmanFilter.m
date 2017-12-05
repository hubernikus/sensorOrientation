function [model, x_filt] ...
            = init_kalmanFilter(modelType, dt_kf, r_circ, omega0, sigma_GPS)
        
model.type = modelType;

switch modelType
    case 'constVel'
        model.dim = 4;

        % A priori Info Conditions
        model.sigma_x = 10;% [m]
        model.sigma_v = 0.10;% [m]

        model.F = [0 0 1 0;   
             0 0 0 1;
             0 0 0 0;
             0 0 0 0];

        model.G = [0 0;
             0 0;
             1 0;
             0 1];

        model.Phi = eye(model.dim) + dt_kf*model.F;

        model.sigma_vt = 0.035; % Uncertainty in Model
        model.a = 0; % uniform constant motion
        model.sigma_at = 0;
        qvn = model.sigma_vt^2; qve = model.sigma_vt^2;
        model.W = diag([qvn, qve]);

        simga_n_GPS = sigma_GPS(1);
        simga_e_GPS = sigma_GPS(2);
        model.R = diag([simga_n_GPS^2, simga_e_GPS^2]);

        % State 
        %x = [r_circ; 0; x_simu(1,1);x_simu(2,1)];
        % x = [x, y, v_x, v_y];
        x_filt = [r_circ; 0; 0; omega0*r_circ];

        model.P0 = diag([model.sigma_x^2, model.sigma_x^2, ...
                         model.sigma_v^2, model.sigma_v^2]);

        model.H = [1 0 0 0;
                   0 1 0 0];
               
    case 'constAcc'
        model.dim = 6;

        % A priori Info Conditions
        model.sigma_p = 10;       % [m]
        model.sigma_v = 0.10;     % [m]
        model.sigma_a = 0.10;     % [m]

        model.F = [0 0 1 0 0 0; 
             0 0 0 1 0 0;
             0 0 0 0 1 0;
             0 0 0 0 0 1;
             0 0 0 0 0 0;
             0 0 0 0 0 0];

        model.G = [0 0;
             0 0;
             0 0;
             0 0;
             1 0;
             0 1];

        model.Phi0 = eye(6) ...
                    + [zeros(4,2), eye(4)*dt_kf; zeros(2,6)] ...
                    + [zeros(2,4), eye(2)*dt_kf^2/2; zeros(4,6)];

        % Uncertainty of motion model
        model.sigma_at = 0.001; 
        model.a_dot = 0; % uniform constant acceleration
        qvn = model.sigma_at^2; qve = model.sigma_at^2;
        model.W = diag([qvn, qve]);
        
        simga_n_GPS = sigma_GPS(1);
        simga_e_GPS = sigma_GPS(2);
        % R [l x l] - Covariance matrix of measurment
        model.R = diag([simga_n_GPS^2, simga_e_GPS^2]);

        % x [n x 1]- State vector
        % x = [x, y, v_x, v_y, a_x, a_y];
        x_filt = [r_circ; 0; 0; omega0*r_circ;  -omega0^2*r_circ; 0];

        % z [l x 1] - Measurements 

        % H [l x n]- Design Matrix ( z = H*x + v )
        model.H = [1 0 0 0 0 0;
                    0 1 0 0 0 0];

        % P [n x n]- Covariance Matrix of the state vector
        model.P0 = diag([model.sigma_p^2, model.sigma_p^2, ...
                        model.sigma_v^2, model.sigma_v^2,...
                        model.sigma_a^2, model.sigma_a^2]);
        
    case 'circularMotion'
        model.dim = 3;

        % Transformation Cartestian - Polar

        % A priori Info Conditions
        model.sigma_p = 10;       % [m]
        model.sigma_v = 0.10;     % [m]
        
        model.F = [0 0 0;
                   0 0 1;
                   0 0 0];

        model.G = [1 0;
                   0 0;
                   0 1];

        model.Phi0 = [1 0 0;
                      0 1 dt_kf;
                      0 0 0];

        
        % Uncertainty of motion model
        model.sigma_rt = 0.001; 
        model.sigma_wt = omega0/20; 
        
        qvn = model.sigma_rt^2; qve = model.sigma_wt^2;
        model.W = diag([qvn, qve]);
        %model.W = model.M^T * model.W  * model.M;
        
        simga_n_GPS = sigma_GPS(1);
        simga_e_GPS = sigma_GPS(2);
        
        % R [l x l] - Covariance matrix of measurment
        model.R = diag([simga_n_GPS^2, simga_e_GPS^2]);

        % x [n x 1]- State vector
        % x = [r, psi, \dot psi];
        x_filt = [r_circ, omega0, 0]';

        % z [l x 1] - Measurements 

        % H [l x n]- Design Matrix ( z = H*x + v )
        model.H = [1 0 0;
                   0 1 0];
               
        % P [n x n]- Covariance Matrix of the state vector
        model.P0 = diag([10^2 (pi/100)^2 (pi/1000)^2]);
        
    otherwise
        fprintf('Not implemented \')
end
        
end