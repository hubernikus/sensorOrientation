function [ DATA ] = gaussMarkovGen( N, corrTime )
% GAUSSMARKOVGENERATOR
% Gauss-Markov process of 1st order x_(k+1) = exp(-1/T_)

whiteNoise = whiteNoiseGen(N);

DATA = zeros(size(whiteNoise));
DATA(1,:) = whiteNoise(1,:);

for jj = 2:N
    DATA(jj,:) = exp(-1/corrTime)*DATA(jj-1,:) + whiteNoise(jj,:);
end

end

