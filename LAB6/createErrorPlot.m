function [max_posErrors, max_velErrors, max_angErrors] = ...
    createErrorPlot(x_sim, x_real, v_sim, v_real, phi_sim, phi_real, titleName)

phi_real = phi_real*180/pi;

err_pos = x_sim-x_real;
absErr_pos = sqrt(sum(err_pos.^2,1));
max_posErrors = max(absErr_pos);

err_vel = v_sim(:,1:end-1)-v_real;
absErr_vel = sqrt(sum(err_vel.^2,1));
max_velErrors = max(absErr_vel);

err_phi = phi_sim*180/pi-phi_real;
absErr_azi = abs(err_phi);
max_angErrors = max(absErr_azi);

if nargin>6
    figure('Position',[0,0,800,1200]); % Plot results
    set(groot,'DefaultAxesFontSize',14)
    set(groot,'DefaultLineLineWidth',1.2)
    
    h1 = subplot(3,1,1);
    plot(h1, phi_real(1:end), absErr_pos); hold on;
    plot(h1, phi_real(1:end), err_pos(1,:)); 
    plot(h1, phi_real(1:end), err_pos(2,:)); 
    legend('Absolut','East direction', 'North direction','Location','northwest')
    ylabel('Position error [m]')
    xlim([phi_real(1),phi_real(end)]);
    
    h2 = subplot(3,1,2);
    plot(h2, phi_real(1:end-1), absErr_vel); hold on;
    plot(h2, phi_real(1:end-1), err_vel(1,:)); 
    plot(h2, phi_real(1:end-1), err_vel(2,:)); 
    legend('Absolut','East direction', 'North direction','Location','northwest')
    ylabel('Velocity error [m]')
    xlim([phi_real(1),phi_real(end-1)]); hold off;
    
    h3 = subplot(3,1,3);
    plot(phi_real, absErr_azi); hold on;
    ylabel('Azimuth Error [deg]'); xlabel('Azimuth [deg]')
    xlim([phi_real(1),phi_real(end)])

    print(strcat('fig/',titleName),'-depsc')
end
