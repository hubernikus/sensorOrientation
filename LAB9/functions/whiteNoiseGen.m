function [ DATA ] = whiteNoiseGen(N, M, var)
% WHITENOSIEGENERATOR 
% White noise: x_(k+1) = w_k

if nargin<2
    M =1;
    if nargin<3
        var = 1;
    end
end

DATA = randn(M,N)*var;

end

