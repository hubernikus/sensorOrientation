%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Sensor Orientation - EPFL 
%
%           Author: Huber Lukas
%           Date: 2017-11-03
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Blackbaord - Comments
%
% WN_ \sigma_{\omega_{WN}} ., \simga_{a_{WN}}
% RC \rightarrow chosen arbitrary \pm 3 \sigma_{RC}
% GM1:  \sigma_{GM}, \beta = \frac 1 T
% b_{k+1} = \Sigma b_k ' + w_{k+1} \quad \text{with} \quad 0 < \Sigma < 1
%
% \sigma_{GM}^2 = \SIGMA ^2 \sigma_{GM}^2 + q
% (1-e^{-12\beta \Delta t}) \sigma_{GM}^2 = q
% \sqrt{q} = \simga_{WN}
% 
%
% Gyro:
% s^{obs} = s^{how} + b_{RC} + b_{GM} + n_{WN}
%
% Accelerometer
% f^b = f^b_{non} + 
% \begin{bmatrix} f_1 RC \\ f_2 RC \end{bmatrix} +
% \begin{bmatrix} randn\\ randn \end{bmatrix} \cdot \sigma_a^2 
%
%
% deg/rn  \qquad 0.05 \deg/\sqrt(n) \quad \rightarrow \quad rad/s/sample
% \frac \pi{180} \cdot \frac 1 \sqrt{3600} \quad \rightarrow \quad rad/ \sqrt(s) = rad/s/ \sqrt{Hz}
% * \sqrt(f_s} \quad \rightarrow \quad rad/s/sam
% 
%%

close all; clear all; clc;

addpath(genpath('../LAB1'))
addpath(genpath('../LAB3'))
%

%% Rescaling sensor error
g = 9.81; % [m/s^2] - Gravity

sampilingFreq = 10;     % [Hz]
%sampilingFreq = 100;    % [Hz]

b_g = 10,               % [deg/h]
simga_G_GM = 0.005;     % [deg/s/sqrt(Hz)]
beta_G_inv = 100;       % [s]
sigma_G_vel = 0.1;      % [deg/sqrt(h)]
b_A = 1;                %[mg]
sigma_A_vel = 50;       % [mu g/ /sqrt(Hz)]

% Rescale to SI units
b_g = b_g/180*pi/3600 
simga_G_GM = simga_G_GM/180*pi
beta_G_inv = beta_G_inv
sigma_G_vel = sigma_G_vel/180*pi*sqrt(samplingFreq)
b_A = b_A * g 
sigma_A_vel = sigma_A_vel / 10e6 * g *sqrt(samplingFreq)
%%

sampleRate = 100; % [Hz]
r_circ = 500; % [m]
omega0 = pi/100; % Frequency rad/s

phi0 = pi/100; % 
azim_init = pi/100; % 


% Simulate Measurement
[acc, gyro, x_real, v_real, phi_real, time] = ... 
           IMUsens_simulation(sampleRate, omega0, r_circ, phi0, azim_init);


% Noise Functions
% whiteNoiseGen()
% randWalkGen()
% gaussMarkovGen()