function [Pxx, freq] = powerSpectralDensity(x, Fs)
%PSD - POWER SPECTRAL DENISTY
%   Caclculate Power spectral Density of the signal x with frequency Fs
%
% -------------------------------------------------------------------------
% Input:
%  x  [Nx1] - Input signal
%  Fs - Frequence of signal (Default: Fs = 1)
%
% -------------------------------------------------------------------------
% Output:
%  freq [N/2+1 x 1] - Frequence values 
%  psdx [N/2+1 x 1] - Power spectral density at correlsponding values 
%
% -------------------------------------------------------------------------

if nargin < 2
    Fs = 1; % Default frequency of 1 Hz
end

% Length of signal
N = length(x);

Pxx = pwelch(x,hanning(N*0.4));

% Correpsonding frequencies
freq = [Fs/length(Pxx):Fs/length(Pxx):Fs]';

end

