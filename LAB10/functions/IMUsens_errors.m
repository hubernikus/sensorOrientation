function [acc_simu, gyro_simu, err] = ...
            IMUsens_errors(acc, gyro, samplingFrequency, g)

g = 9.81; %Gravitiy [m/s^2]

N_samples = length(gyro);

err  = [];
%% Initialization
% Gyroscope error
gyro_bias = 18; %deg/h
gyro_GM1.sigma = 0.01; % [rad/s/sqrt(Hz)]
gyro_GM1.invBeta = 30; % [s]
gyro_rw = 0.1; % [deg/sqrt(h)]


err.gyro_bias = gyro_bias /180*pi; % [rad/s]
err.gyro_GM1.initSigma = 0;
err.gyro_GM1.sigma = gyro_GM1.sigma /180*pi * sqrt(samplingFrequency); % [rad/s/sqrt(Hz)]
err.gyro_GM1.beta = 1/gyro_GM1.invBeta; % [Hz]
err.gyro_rw = gyro_rw /180*pi *1/sqrt(3600) *sqrt(samplingFrequency); % [rad/sqrt(s)]



% Accelerometer error
acc_wn = 50;  % [mu*g/sqrt(Hz]
%acc1_GM1 =[];
acc1_GM1.initSigma= -100; % [mu*g]
acc1_GM1.sigma = 50; % [mu*g/sqrt(Hz)]
acc1_GM1.invBeta = 60; % [s]

%acc2_GM1 =[];
acc2_GM1.initSigma = 200; % [mu*g]
acc2_GM1.sigma = 50; % [mu*g/sqrt(Hz)]
acc2_GM1.invBeta = 60; % [s]


err.acc_wn = acc_wn*1e-6*g*sqrt(samplingFrequency);  % [m/s^2]
err.acc1_GM1.initSigma= acc1_GM1.initSigma*1e-6*g; % [m/s^2]
err.acc1_GM1.sigma = acc1_GM1.sigma*1e-6*g*sqrt(samplingFrequency); % [m/s^2/qr(Hz)]
err.acc1_GM1.beta = 1/acc1_GM1.invBeta; % [Hz]

err.acc2_GM1.initSigma = acc2_GM1.initSigma*1e-6*g; % [m/s^2]
err.acc2_GM1.sigma = acc2_GM1.sigma*1e-6*g*sqrt(samplingFrequency); % [m/s^2/qr(Hz)]
err.acc2_GM1.beta = 1/acc2_GM1.invBeta; % [Hz]

%% Apply Errors
dim = 1;

% Gyro Bias
gyro_bias = ones(dim,N_samples)*err.gyro_bias;

% Gyro gauss markov
gyro_GM =  gaussMarkovGen(N_samples, dim, samplingFrequency, err.gyro_GM1.beta, err.gyro_GM1.sigma);

% Gyro random walk
gyro_randomWalk =  randWalkGen(N_samples, dim, err.gyro_rw);


gyro_simu = gyro + gyro_bias + gyro_GM + gyro_randomWalk;


% Accelerometer white noise
dim =2;
acc_noise = whiteNoiseGen(N_samples, dim, err.acc_wn);

% Accelerometer Gauss Markov
dim = 1;
acc1_GM =  gaussMarkovGen(N_samples, dim, samplingFrequency, err.acc1_GM1.beta, err.acc1_GM1.sigma, err.acc1_GM1.initSigma);
acc2_GM =  gaussMarkovGen(N_samples, dim, samplingFrequency, err.acc2_GM1.beta, err.acc2_GM1.sigma, err.acc2_GM1.initSigma);

acc_simu = acc + acc_noise + [acc1_GM;acc2_GM];

end