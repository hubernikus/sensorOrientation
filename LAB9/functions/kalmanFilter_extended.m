function [x_filt, innovation, sigma_pred, x_tild] ...
            = kalmanFilter_extended(model, x_simu, x_filt,dt_kf, dt_gps)
% Extract variables
P = model.P0;
F = model.F;
G = model.G;
dim = model.dim;
R = model.R;
W = model.W;
H = model.H;

% Initialization
% Auxiliary matrix A
A = [-F, G*W*G'; ...
     zeros(dim), F'] * dt_kf;

B = expm(A);

Phi = B(dim+1:2*dim,dim+1:2*dim)';
Q_k = Phi*B(1:dim,dim+1:2*dim);

% Prediction
x_tild(:,1) = Phi*x_filt(:,1);
P_tilde = Phi*P*Phi' + Q_k;

N_sim = size(x_simu,2);

for ii = 1:N_sim
    % Mesurement --- z, R ----

    % Gain (weight)
    K = P_tilde*H'/(H*P_tilde*H' + R);

    % State update
    innovation(:,ii) = (x_simu(:,ii)- H*x_tild(:,1+((ii-1)*dt_gps/dt_kf) ));
    x_filt(:,ii+1) = x_tild(:,1+((ii-1)*dt_gps/dt_kf)) + K*innovation(:,ii);
    
    
    % Covariance update
    P = (eye(dim)- K*H)*P_tilde;

    % KF-predicted
    sigma_pred(ii) = sqrt(sum(diag(P)));

    % Prediction - Inner Loop
    % Output ---- x_est, P ---

    % Auxiliary matrix A

    for jj = 1:dt_gps/dt_kf    

        % Output ---- x_est, P ---

        % Auxiliary matrix A
        A = [-F, G*W*G'; ...
             zeros(dim), F'] * dt_kf;

        B = expm(A);

        Phi = B(dim+1:2*dim,dim+1:2*dim)';
        Q_k = Phi*B(1:dim,dim+1:2*dim);

        % Prediction
        if(jj==1) % first round take measurement
            x_tild(:,jj+1 + (ii-1)*dt_gps/dt_kf) = Phi*x_filt(:,ii+1);
            P_tilde = Phi*P*Phi' + Q_k;
        else    % take prediction after
            x_tild(:,jj+1 + (ii-1)*dt_gps/dt_kf) = Phi*x_tild(:,jj+(ii-1)*dt_gps/dt_kf);
            P_tilde = Phi*P_tilde*Phi' + Q_k;
        end
    end
end
end