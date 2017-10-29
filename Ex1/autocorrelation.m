function [phi, xAxis, sigma,T] = autocorrelation(z)
% AUTOCORELLATION - CACLULATION OF THE AUTOCORRELATION OF A FUNCION
%
% -------------------------------------------------------------------------
% INPUT:
%  z - input signal
%
% -------------------------------------------------------------------------
% OUTPUT:
%  phi - Autocorrelation 
%
% -------------------------------------------------------------------------

N = length(z); % Legnth of the vector
m = 1/N*sum(z); % Mean value of z

phi = xcorr(z,'unbiased');

sigma = sqrt(phi(N));

phi = phi/phi(N); % Normalize

xAxis = (-(N-1):N-1)'; % Normalize

% Find the position where the autocorrelation passed the value e^-1
pass_e1 = find(phi(N+1:end)< exp(-1),1); % symetric, first value on the left is found.
T = pass_e1

end