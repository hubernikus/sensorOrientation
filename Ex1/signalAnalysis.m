function [ sigma_corr, sigma_sig ,pass_e1] = signalAnalysis( DATA , figName )

% Get data dimensions
[M,N] = size(DATA);

% De-mean data
DATA = DATA - mean(DATA,1);

%% Generation of 3 random white noice sequences
cols = ['r','g','b','y','c'];

%
if(N>length(cols)); warning('Not enough color specifications!'); end

figure('Position',[0,0,800,1200]); subplot(3,1,1)  % Plot results
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.2)
for i=1:N
    plot(DATA(:,i),cols(i)); hold on;
end
xlabel('Sample Number []'); ylabel('Signal value []')
xlim([0,M])


subplot(3,1,2)
sigma_corr = zeros(3,1);
pass_e1= zeros(3,1);

for i = 1:N
    % Autocorrelation Using User defined Formulla
    [phi, xAxis, sigma_corr(i),pass_e1] = autocorrelation(DATA(:,i));
    plot(xAxis,phi, cols(i)); hold on;   
end
plot([xAxis(1), xAxis(end)],[exp(-1),exp(-1)],'k--')
xlim([xAxis(1), xAxis(end)]);ylim([-1.2,1.2])
xlabel('Sample Number []'); ylabel('Autocorrelation []')

subplot(3,1,3)
yRange = [1,1];
Fs = 1; % Default frequency: 1Hz
for i = 1:N
    % PSD Using User defined Formulla
    [psd_DATA, freq1] = powerSpectralDensity(DATA(:,i), Fs);
    loglog(freq1, psd_DATA, cols(i)); hold on;
    yRange(1) = min(yRange(1), min(psd_DATA));
    yRange(2) = max(yRange(2), max(psd_DATA));
end
xlabel('Frequency [log]'); ylabel('Power Spectral Density [log]')
xlim([min(freq1),max(freq1)]); 
ylim([yRange(1)*0.9, yRange(2)*1.1]);


figure;
subplot(4,1,4)
yRange = [1,1];
Fs = 1; % Default frequency: 1Hz
for i = 1:N
    % PSD Using User defined Formulla
    allan_DATA= allandev(DATA(:,i),'AllanDev');
    loglog(allan_DATA, cols(i)); hold on;
    yRange(1) = min(yRange(1), min(allan_DATA));
    yRange(2) = max(yRange(2), max(allan_DATA));
end
xlabel('Frequency [log]'); ylabel('Power Spectral Density [log]')
xlim([min(freq1),max(freq1)]); 
ylim([yRange(1)*0.9, yRange(2)*1.1]);



%print(strcat('fig/',figName),'-dpdf')
print(strcat('fig/',figName),'-depsc')

sigma_sig = std(DATA); %% change if wanted




end

