function [ sigma ] = signalAnalysis( DATA , figName )

%% Generation of 3 random white noice sequences
cols = ['r','g','b'];

figure('Position',[0,0,800,1200]); subplot(3,1,1)  % Plot results
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)
for i=1:3
    plot(DATA(:,i),cols(i)); hold on;
end
xlabel('Sample Number []'); ylabel('Signal value []')


subplot(3,1,2)
sigma = zeros(3,1);
for i = 1:3
    % Autocorrelation Using User defined Formulla
    [phi, xAxis, sigma(i)] = autocorrelation(DATA(:,i));
    Nc = length(phi);
    plot(xAxis,phi, cols(i)); hold on;   
end
plot([xAxis(1), xAxis(end)],[exp(-1),exp(-1)],'k--')
ylim([-1.2,1.2])
xlabel('Sample Number []'); ylabel('Autocorrelation []')

subplot(3,1,3)
yRange = [1,1];
Fs = 1; % Default frequency: 1Hz
for i = 1:3
    % PSD Using User defined Formulla
    [psd_DATA, freq1] = powerSpectralDensity(DATA(:,i), Fs);
    loglog(freq1, psd_DATA, cols(i)); hold on;
    yRange(1) = min(yRange(1), min(psd_DATA));
    yRange(2) = max(yRange(2), max(psd_DATA));
end
xlabel('Frequency [log]'); ylabel('Power Spectral Density [log]')
xlim([min(freq1),max(freq1)]); 
ylim([yRange(1)*0.9, yRange(2)*1.1]);

print(strcat('fig/',figName),'-dpdf')
print(strcat('fig/',figName),'-depsc')

sigma = std(DATA); %% change if wanted

end

