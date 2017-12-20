function [ DATA ] = gaussMarkovGen(N_samples, dim, samplingFreq, beta,  var, varInit)
% GAUSSMARKOVGENERATOR
% Gauss-Markov process of 1st order x_(k+1) = exp(-1/T_)

deltaT = 1/samplingFreq;

if nargin < 4
    beta = 1;
    if nargin < 5
        var = 1;
    end
end

whiteNoiseVar = sqrt((1-exp(-2*deltaT*beta))*var^2);

whiteNoise = randn(dim,N_samples)*whiteNoiseVar;

DATA = zeros(size(whiteNoise));

if nargin<6
    DATA(:,1) = whiteNoise(:,1);
else
    DATA(:,1) = randn(dim,1)*varInit;
end

for jj = 2:N_samples
    DATA(:,jj) = exp(-deltaT*beta)*DATA(:,jj-1) + whiteNoise(:,jj);
end

end

