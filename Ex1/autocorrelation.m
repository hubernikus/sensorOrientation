function [phi, xAxis, sigma] = autocorrelation(z)
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
end