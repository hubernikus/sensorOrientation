function [ DATA ] = gaussMarkovGen(N, M, deltaT, corrTime,  var)
% GAUSSMARKOVGENERATOR
% Gauss-Markov process of 1st order x_(k+1) = exp(-1/T_)

if nargin < 3
    corrTime = 1;
    if nargin < 4
        var = 1;
    end
end

whiteNoiseVar = sqrt((1-exp(-2*deltaT/corrTime))*var^2),

whiteNoise = whiteNoiseGen(N, M, whiteNoiseVar);

DATA = zeros(size(whiteNoise));
DATA(:,1) = whiteNoise(:,1);

for jj = 2:N
    DATA(:,jj) = exp(-deltaT/corrTime)*DATA(:,jj-1) + whiteNoise(:,jj);
end

end

