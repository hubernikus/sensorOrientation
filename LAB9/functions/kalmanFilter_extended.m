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

M = eye(2);

% Initialization
% Auxiliary matrix A
if strcmp(model.type,'circularMotion') % transform radiational - cartesian
    x_filt_cart = convPolCart(x_filt);
    
    M = mCalc(x_filt_cart(1),x_filt_cart(2));
    W = M*W*M';
    
end

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
    innovation(:,ii) = (M*x_simu(:,ii)- H*x_tild(:,1+((ii-1)*dt_gps/dt_kf) ));
    x_filt(:,ii+1) = x_tild(:,1+((ii-1)*dt_gps/dt_kf)) + K*innovation(:,ii);
    if strcmp(model.type,'circularMotion') % transform radiational - cartesian
        x_filt_cart(:,ii+1) = convPolCart(x_filt(:,ii+1));

        M = mCalc(x_filt_cart(1,ii+1),x_filt_cart(2,ii+1));
        W = M*W*M';
    end
    
    
    
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
        if strcmp(model.type,'circularMotion') % transform radiational - cartesian
            M = mCalc(x_filt_cart(1,ii+1),x_filt_cart(2,ii+1));
            W = M*W*M';
        end
        
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

if strcmp(model.type,'circularMotion') % transform radiational - cartesian
    x_filt = x_filt_cart; % return x_filt in cartesian coordinates
end

end

function M = mCalc(p_n, p_e)
r = sqrt(p_n^2 + p_e^2);
psi = atan2(p_e,p_n);
M = [cos(psi), sin(psi);
    -1/r*sin(psi), 1/r*cos(psi)];
end

function M = mCalcAng(psi)
M = [cos(psi), sin(psi);
    -1/r*sin(psi), 1/r*cos(psi)];
end

function x_cart = convPolCart(x_pol)
x_cart = [x_pol(1)*cos(x_pol(2)); % p_n
          x_pol(1)*sin(x_pol(2)); % p_e
          -x_pol(1)*x_pol(3)*sin(x_pol(2));  % v_n
          -x_pol(1)*x_pol(3)*cos(x_pol(2))]; % v_e
      
end

function x_pol = convCartPol(x_cart)
x_pol = [sqrt(x_cart(1)^2 + x_cart(2)^2)
         atan2(x_cart(2),x_cart(1));
         0];
end