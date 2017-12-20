function [model, x_filt] ...
            = standardDeviation_pos(modelType, dt_kf, r_circ, omega0, sigma_GPS)
        
iter = 1:size(x_filt,3)        
    x_filt = x_KF(:,:,iter);
    % Empirical standard deviaion characterizing GPS
    sigma_gps_x = std(x_simu(1,:)-x_real(1,:));
    sigma_gps_y = std(x_simu(2,:)-x_real(2,:));
    sigma_gps(iter) = sqrt(sigma_gps_x^2+sigma_gps_y^2);

    % Empirical standard deviaion characterizing GPS
    sigma_filt_x = std(x_filt(1,1:end-1)-x_real(1,:));
    sigma_filt_y = std(x_filt(2,1:end-1)-x_real(2,:));
    sigma_filt(iter) = sqrt(sigma_filt_x^2+sigma_filt_y^2);

    % Empirical standard deviaion characterizing GPS
    meanSqr_gps(iter) = sum(sigma_gps_x^2+sigma_gps_y^2);

    % Empirical standard deviaion characterizing GPS
    meanSqr_filt(iter) = sum(sigma_filt_x^2+sigma_filt_y^2);

end