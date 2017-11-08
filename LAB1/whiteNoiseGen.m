function [ DATA ] = whiteNoiseGen( N, M, var )
% WHITENOSIEGENERATOR 
% White noise: x_(k+1) = w_k

if nargin<3
    var = 1;
    if nargin<2
        M = 3;
    end
end

DATA = randn(N,1)*var;

end

