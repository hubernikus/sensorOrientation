function [ DATA ] = whiteNoiseGen( N )
% WHITENOSIEGENERATOR 
% Random walk (RW x_(k+1) = x_k + w_k)

whiteNoise = randn(N,3);
DATA = cumsum(whiteNoise,1);

end

