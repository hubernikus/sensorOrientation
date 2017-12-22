function [X, dX] = inertialNavigationStep(X0, acc, gyr, dt)
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

alpha0 = X0(1);
v0 = X0(2:3);
x0 = X0(4:5);

%acc(:) = acc(:) - X0(8:9);
%gyr(:) = gyr(:) - X0(6)-X0(7);

R_b_m = [cos(alpha0), -sin(alpha0); sin(alpha0), cos(alpha0)];
 
switch intOrder
    case 1
        % rapezodial attitude approximation
        alpha = alpha0 + gyr*dt;

        Omega_mb_b = [0, -gyr; gyr, 0];
    
        R_b_m = R_b_m*expm(Omega_mb_b*dt); 
        %R_b_m = [cos(alpha(k)), - sin(alpha(k)); sin(alpha(k)), cos(alpha(k))];

        % Navigation computation
        v_m = v0 + R_b_m*acc*(dt);
        x_m = x0 + v_m*dt;
    case 2
        % rapezodial attitude approximation
        alpha = alpha0 + gyr*dt;

        Omega_mb_b = [0, -gyr; gyr, 0];
    
        R_b_m = R_b_m*expm(Omega_mb_b*dt); 
        %R_b_m = [cos(alpha(k)), - sin(alpha(k)); sin(alpha(k)), cos(alpha(k))];

        % Navigation computation
        v_m = v0 + R_b_m*acc*(dt);
        x_m = x0 + v_m*dt;
        
        R_b_m =  {R_b_m,  R_b_m};
        % rapezodial attitude approximation
        alpha(k) = alpha(k-1) + 0.5*(gyr(k)+gyr(k-1))*(time(k)-time(k-1));

        Omega_mb_b = [0 -gyr(k); gyr(k) 0];

        R_b_m{1} = R_b_m{2};  
        R_b_m{2} = R_b_m{1}*expm(Omega_mb_b*(time(k)-time(k-1))); 
        %R_b_m = [cos(alpha(k)), - sin(alpha(k)); sin(alpha(k)), cos(alpha(k))];

        % Navigation computation
        v_m(:,k) = v_m(:,k-1) + 0.5*(R_b_m{2}*acc(:,k)+R_b_m{1}*acc(:,k-1))*(time(k)-time(k-1));
        x_m(:,k) = x_m(:,k-1) + 0.5*(v_m(:,k)+v_m(:,k-1))*(time(k)-time(k-1));
    otherwise
        warning('Method not yet implemented')
end

X = [alpha; v_m; x_m; X0(6:9)];
dX = X-X0;

end
% 
% if(length(x0) == 2) % 2D case
%     
%     % Initial rotation
%     alpha = alpha0; 
%     
%     % Initial position
%     v_m(:,1) = v0;
%     x_m(:,1) = x0;
%     
%     R_b_m = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
%             
%     switch intOrder
%         case 1
%             % Dead reckoning
%             for k = 2:length(time) % first timestep at 0 - ind 1
%                 % rapezodial attitude approximation
%                 alpha(k) = alpha(k-1) + 1*gyr(k)*(time(k)-time(k-1));
% 
%                 Omega_mb_b = [0 -gyr(k); gyr(k) 0];
% 
%                 R_b_m = R_b_m*expm(Omega_mb_b*(time(k)-time(k-1))); 
%                 %R_b_m = [cos(alpha(k)), - sin(alpha(k)); sin(alpha(k)), cos(alpha(k))];
% 
%                 % Navigation computation
%                 v_m(:,k) = v_m(:,k-1) + R_b_m*acc(:,k)*(time(k)-time(k-1));
%                 x_m(:,k) = x_m(:,k-1) + v_m(:,k)*(time(k)-time(k-1));
%             end
%         
%         case 2
%             R_b_m =  {R_b_m,  R_b_m};
%             % Dead reckoning
%             for k = 2:length(time) % first timestep at 0 - ind 1
%                 % rapezodial attitude approximation
%                 alpha(k) = alpha(k-1) + 0.5*(gyr(k)+gyr(k-1))*(time(k)-time(k-1));
% 
%                 Omega_mb_b = [0 -gyr(k); gyr(k) 0];
% 
%                 R_b_m{1} = R_b_m{2};  
%                 R_b_m{2} = R_b_m{1}*expm(Omega_mb_b*(time(k)-time(k-1))); 
%                 %R_b_m = [cos(alpha(k)), - sin(alpha(k)); sin(alpha(k)), cos(alpha(k))];
% 
%                 % Navigation computation
%                 v_m(:,k) = v_m(:,k-1) + 0.5*(R_b_m{2}*acc(:,k)+R_b_m{1}*acc(:,k-1))*(time(k)-time(k-1));
%                 x_m(:,k) = x_m(:,k-1) + 0.5*(v_m(:,k)+v_m(:,k-1))*(time(k)-time(k-1));
%             end
% 
%         otherwise
%             warning('Method not yet implemented')
%     end
% 
% else
%     error('3D not yet implemented \n')
% end
% 
% phi = alpha-pi;
% 
% end