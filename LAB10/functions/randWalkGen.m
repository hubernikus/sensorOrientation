function [ DATA ] = randWalkGen( N, dim , std)
% WHITENOSIEGENERATOR 
% Random walk (RW x_(k+1) = x_k + w_k)

if nargin<2
    dim = 1;
    if nargin<3
        std = 1;
    end
end

whiteNoise = randn(dim, N)*std;
DATA = cumsum(whiteNoise,1);

end