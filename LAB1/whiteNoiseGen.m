function [ DATA ] = whiteNoiseGen( N, M var )
% WHITENOSIEGENERATOR 
% White noise: x_(k+1) = w_k

if nargin<3
    M =1;
    if nargin<2
        var = 1;
    end
end

DATA = randn(N,M)*var;

end

