function [x_m, v_m, phi] = inertialNavigation(x0, v0, alpha0, acc, gyr, time, intOrder)
% Straspdown intertial navigation0
%
% x0 - initial position
% v0 - initial velocity
% alph0 - initial orientation0
%
% acc - accelerometer measurments
% gyr - gyroscope measurments
%
% time - simulation time [s]
%
% intOrder - Integration order 1 2 6

if nargin<7
    intOrder = 1;
end

if(length(x0) == 2) % 2D case
    
    % Initial rotation
    alpha = alpha0; 
    
    % Initial position
    v_m(:,1) = v0;
    x_m(:,1) = x0;
    
    R_b_m = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
            
    switch intOrder
        case 1
            % Dead reckoning
            for k = 2:length(time) % first timestep at 0 - ind 1
                % rapezodial attitude approximation
                alpha(k) = alpha(k-1) + 1*gyr(k)*(time(k)-time(k-1));

                Omega_mb_b = [0 -gyr(k); gyr(k) 0];

                R_b_m = R_b_m*expm(Omega_mb_b*(time(k)-time(k-1))); 
                %R_b_m = [cos(alpha(k)), - sin(alpha(k)); sin(alpha(k)), cos(alpha(k))];

                % Navigation computation
                v_m(:,k) = v_m(:,k-1) + R_b_m*acc(:,k)*(time(k)-time(k-1));
                x_m(:,k) = x_m(:,k-1) + v_m(:,k)*(time(k)-time(k-1));
            end
        
        case 2
            R_b_m =  {R_b_m,  R_b_m};
            % Dead reckoning
            for k = 2:length(time) % first timestep at 0 - ind 1
                % rapezodial attitude approximation
                alpha(k) = alpha(k-1) + 0.5*(gyr(k)+gyr(k-1))*(time(k)-time(k-1));

                Omega_mb_b = [0 -gyr(k); gyr(k) 0];

                R_b_m{1} = R_b_m{2};  
                R_b_m{2} = R_b_m{1}*expm(Omega_mb_b*(time(k)-time(k-1))); 
                %R_b_m = [cos(alpha(k)), - sin(alpha(k)); sin(alpha(k)), cos(alpha(k))];

                % Navigation computation
                v_m(:,k) = v_m(:,k-1) + 0.5*(R_b_m{2}*acc(:,k)+R_b_m{1}*acc(:,k-1))*(time(k)-time(k-1));
                x_m(:,k) = x_m(:,k-1) + 0.5*(v_m(:,k)+v_m(:,k-1))*(time(k)-time(k-1));
            end

        otherwise
            warning('Method not yet implemented')
    end

else
    error('3D not yet implemented \n')
end

phi = alpha-pi;

end