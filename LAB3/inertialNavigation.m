function [x_m, v_m, alpha] = inertialNavigation(x0, v0, alpha0, acc, gyr, time, intOrder)
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

            % Atttitude Computation
            %dR_b_m = R_b_m*Omega_mb_b;
            
            % Assumption of constant angular rate
            
            %R_b_m = R_b_m*expm(Omega_mb_b*(time(k)-time(k-1))) % HOW TO!!!
            R_b_m = [cos(alpha(k)), - sin(alpha(k)); sin(alpha(k)), cos(alpha(k))];

            % Navigation computation
            v_m(:,k) = v_m(:,k-1) + R_b_m*acc(:,k)*(time(k)-time(k-1));
            x_m(:,k) = x_m(:,k-1) + v_m(:,k)*(time(k)-time(k-1));
            
        end
        
        case 2
        % Dead reckoning
        for k = 2:length(time) % first timestep at 0 - ind 1
            % rapezodial attitude approximation
            alpha(k) = alpha(k-1) + 0.5*(gyr(k)+gyr(k-1))*(time(k)-time(k-1));
            
            Omega_mb_b = [0 -gyr(k); gyr(k) 0];

            % Atttitude Computation
            %dR_b_m = R_b_m*Omega_mb_b;
            
            % Assumption of constant angular rate
            
            %R_b_m = R_b_m*expm(Omega_mb_b*(time(k)-time(k-1))) % HOW TO!!!
            R_b_m = [cos(alpha(k)), - sin(alpha(k)); sin(alpha(k)), cos(alpha(k))];

            % Navigation computation
            v_m(:,k) = v_m(:,k-1) + R_b_m*0.5*(acc(:,k)+acc(:,k-1))*(time(k)-time(k-1));
            x_m(:,k) = x_m(:,k-1) + 0.5*(v_m(:,k)+v_m(:,k-1))*(time(k)-time(k-1));
        end
%         case 2
%         
%         i = 2;
%         % Dead reckoning
%         alpha(:,i) = alpha(:,i-1) + gyr(i)*(time(i)-time(i-1));
%         %alpha = alpha + 0.5*(gyr(i+1)-gyr(i))*(dT);
% 
%         
%         R_b_m = [cos(alpha(:,i)), -sin(alpha(:,i)); sin(alpha(:,i)), cos(alpha(:,i))];
%         Omega_mb_b = [0 -gyr(i); gyr(i) 0];
% 
%         % Atttitude Computation
%         %dR_b_m = R_b_m(:,:,1)*Omega_mb_b;
% 
%         % Assumption of constant angular rate
%         R_b_m(:,:,2) = exp(Omega_mb_b*(time(i)-time(i-1)));
% 
%         % Navigation computation
%         v_m(:,i) = v_m(:,i-1) + R_b_m*acc(:,i-1)*(time(i)-time(i-1));
%         x_m(:,i) = x_m(:,i-1) + v_m(:,i)*(time(i)-time(i-1));
%         
%         for i = 3:length(time) % first timestep at 0 - ind 1
%             % rapezodial attitude approximation
%             alpha(:,i) = alpha(:,i-1) + 0.5*(gyr(i)+gyr(i-1))*(time(i)-time(i-1));
%             
%             %R_b_m(:,:,2) = [cos(alpha(:,i)), -sin(alpha(:,i)); sin(alpha(:,i)), cos(alpha(:,i))];
%             Omega_mb_b = [0 -alpha  (i); gyr(i) 0];
% 
%             % Atttitude Computation
%             %dR_b_m = R_b_m*Omega_mb_b;
%             
%             % Assumption of constant angular rate
%             R_b_m = R_b_m*exp(Omega_mb_b*(time(i)-time(i-1)));
% 
%             % Navigation computation
%             v_m(:,i) = v_m(:,i-1) + R_b_m*acc(:,i-1)*(time(i)-time(i-1));
%             x_m(:,i) = x_m(:,i-1) + v_m(:,i)*(time(i)-time(i-1));
%             
%         end
        otherwise
            warning('Method not yet implemented \n')
    end

else
    error('3D not yet implemented \n')
end



end