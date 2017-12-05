function [model, x_filt] ...
            = init_kalmanFilter(modelType, dt_kf, r_circ, omega0, sigma_GPS)
        
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

        model.Phi = eye(dim) + dt_kf*F;

        model.sigma_vt = 0.035; % Uncertainty in Model
        model.a = 0; % uniform constant motion
        model.sigma_at = 0;
        model.qvn = sigma_vt^2; qve = sigma_vt^2;
        model.W = diag([qvn, qve]);

        model.R = diag([simga_x_GPS^2, simga_y_GPS^2]);

        % State 
        %x = [r_circ; 0; x_simu(1,1);x_simu(2,1)];
        % x = [x, y, v_x, v_y];
        x_filt = [0; r_circ; -omega0*r_circ; 0];

        model.P0 = diag([sigma_x^2, sigma_x^2, sigma_v^2, sigma_v^2]);

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
        
        simga_x_GPS = sigma_GPS(1);
        simga_y_GPS = sigma_GPS(2);
        % R [l x l] - Covariance matrix of measurment
        model.R = diag([simga_x_GPS^2, simga_y_GPS^2]);

        % x [n x 1]- State vector
        % x = [x, y, v_x, v_y, a_x, a_y];
        x_filt = [0; r_circ; -omega0*r_circ; 0; 0; -omega0^2*r_circ];

        % z [l x 1] - Measurements 

        % H [l x n]- Design Matrix ( z = H*x + v )
        model.H = [1 0 0 0 0 0;
                0 1 0 0 0 0];

        % P [n x n]- Covariance Matrix of the state vector
        model.P0 = diag([model.sigma_p^2, model.sigma_p^2, ...
                        model.sigma_v^2, model.sigma_v^2,...
                        model.sigma_a^2, model.sigma_a^2]);
        
    case 'realModel'
        fprintf('Not implemented \')
    otherwise
        fprintf('Not implemented \')
end
        
end