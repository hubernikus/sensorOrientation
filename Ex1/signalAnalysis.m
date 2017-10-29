function [ sigma_corr, sigma_sig ,pass_e1] = signalAnalysis( DATA , figName, dt)


if nargin < 3
    dt =1;
end

% Get data dimensions
[M,N] = size(DATA);

% De-mean data
DATA = DATA - mean(DATA,1);

%% Generation of 3 random white noice sequences
cols = ['b','r','g','y','c'];

%
if(N>length(cols)); warning('Not enough color specifications!'); end

figure('Position',[0,0,800,1200]); subplot(3,1,1)  % Plot results
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
for i=1:N
    xDat = (0:(M-1))*dt;
    plot(xDat, DATA(:,i),cols(i)); hold on;
end
xlabel('Time [s]'); ylabel('Signal value []')
xlim([0,xDat(end)])


subplot(3,1,2)
sigma_corr = zeros(N,1);
pass_e1= zeros(N,1);
for i = 1:N
    % Autocorrelation Using User defined Formulla
    [phi, xAxis, sigma_corr(i),pass_e1(i)] = autocorrelation(DATA(:,i));
    xAxis = xAxis*dt;
    plot(xAxis,phi, cols(i)); hold on;   
end
plot([xAxis(1), xAxis(end)],[exp(-1),exp(-1)],'k--')
xlim([xAxis(1), xAxis(end)]);ylim([-1.2,1.2])
xlabel('Time [s]'); ylabel('Autocorrelation []')



subplot(3,1,3)
yRange = [1,1];
Fs = 1; % Default frequency: 1Hz
for i = 1:N
    % PSD Using User defined Formulla
    [psd_DATA, freq1] = powerSpectralDensity(DATA(:,i), Fs);
    freq1 = freq1/dt;
    plot(freq1, psd_DATA, cols(i)); hold on;
    yRange(1) = min(yRange(1), min(psd_DATA));
    yRange(2) = max(yRange(2), max(psd_DATA));
end
%xlabel('Frequency [log]'); ylabel('Power Spectral Density [log]')
xlabel('Frequency [Hz]'); ylabel('Power Spectral Density')
xlim([min(freq1),max(freq1)]); 
%ylim([yRange(1)*0.9, yRange(2)*1.1]);
%ylim([yRange(1)-1, yRange(2)+1]);

print(strcat('fig/',figName),'-depsc')

%print(strcat('fig/',figName),'-dpdf')

sigma_sig = std(DATA)'; %% change if wanted

end

