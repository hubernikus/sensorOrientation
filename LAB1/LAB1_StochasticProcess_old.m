%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      EPFL - Sensor Orientation
%                     LAB 1 - Stochastic Processes
%                              Lukas Huber
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear variables;

% Initialize (pseudo-)random number generator based on birthday
rng(19920526);  

fprintf('Program started. \n');

%% Part A
close all; 
N_random = 200000;  % Length of sequences

%% Generation of 3 random white noice sequences
cols = ['b','g','r'];

DATA = whiteNoiseGen(N_random);

figure('Position',[0,0,800,1200]); subplot(4,1,1)  % Plot results
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)
for i=1:3
    plot(DATA(:,i),cols(i)); hold on;
end
xlabel('Sample Number []'); ylabel('Signal value []')


subplot(4,1,2)
for i = 1:3
    % Autocorrelation Using User defined Formulla
    %phi = autocorrelation(whiteNoise(:,i));
    phi = autocorrelation(DATA(:,i),'unbiased');
    %phi = [flip(phi),phi(2:end)];
    Nc = length(phi);
    plot((1:Nc)-1,phi, cols(i)); hold on;
   
end
ylim([-1.2,1.2])
xlabel('Sample Number []'); ylabel('Autocorrelation []')

subplot(4,1,3)
yRange = [1,1];

Fs = 1; % Default frequency: 1Hz
for i = 1:3
    % PSD Using User defined Formulla
    [psd_randomWalk1, freq1] = powerSpectralDensity(DATA(:,i), Fs);
    loglog(freq1, psd_randomWalk1, cols(i)); hold on;
    yRange(1) = min(yRange(1), min(psd_randomWalk1));
    yRange(2) = max(yRange(2), max(psd_randomWalk1));
end
xlabel('Sample Number [log]'); ylabel('Power Spectral Density [log]')
xlim([min(freq1),max(freq1)]); ylim([yRange(1)*0.9, yRange(2)*1.1]);


% subplot(4,1,4)
% yRange = [1,1];
% Fs = 1; % Default frequency: 1Hz
% for i = 1:3
%     % PSD Using User defined Formulla
%     allanDev_whitNoise = allandev(whiteNoise(:,i), 'whiteNoise');
%     loglog(freq1, psd_randomWalk1, cols(i)); hold on;
%     yRange(1) = min(yRange(1), min(psd_randomWalk1));
%     yRange(2) = max(yRange(2), max(psd_randomWalk1));
% end
% xlabel('Sample Number [log]'); ylabel('Power Spectral Density [log]')
% xlim([min(freq1),max(freq1)]); ylim([yRange(1)*0.9, yRange(2)*1.1]);

%print('fig/WhiteNoise','-dpdf')
print('fig/WhiteNoise','-depsc')


randomWalk = cumsum(DATA, 1);

figure('Position',[0,0,800,1200]); subplot(4,1,1)  % Plot results
set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)
for i=1:3
    plot(randomWalk(:,i),cols(i)); hold on;
end
xlabel('Sample Number []'); ylabel('Signal value []')

subplot(4,1,2); 
yRange = [0,0];
for i = 1:3
    % Autocorrelation Using User defined Formulla
    phi = xcorr(DATA(:,i),'unbiased');
    %phi = autocorrelation(randomWalk(:,i));
    Nc = length(phi);
    plot((1:Nc)-1,phi, cols(i)); hold on;   
    yRange(1) = min(yRange(1), min(phi));
    yRange(2) = max(yRange(2), max(phi));
end
ylim([yRange(1)-10, yRange(2)+10]);
xlabel('Sample Number []'); ylabel('Autocorrelation []')


subplot(4,1,3)
%fonts: set(groot,'DefaultAxesFontSize',14)
%set(groot,'DefaultLineLineWidth',1.5)
yRange = [1,1];

Fs = 1; % Default frequency: 1Hz
for i = 1:3
    % Autocorrelation Using User defined Formulla
    [psd_randWalk1, freq1] = powerSpectralDensity(randomWalk(:,i), Fs);
    loglog(freq1, psd_randWalk1, cols(i)); hold on;
    yRange(1) = min(yRange(1), min(psd_randWalk1));
    yRange(2) = max(yRange(2), max(psd_randWalk1));
end
xlabel('Sample Number [log]'); ylabel('Power Spectral Density [log]')
xlim([min(freq1),max(freq1)]); ylim([yRange(1)*0.9, yRange(2)*1.1]);

print('fig/RandomWalk','-dpdf')
print('fig/RandomWalk','-depsc')


%% Gauss Markov Process  walk (RW x_(k+1) = x_k + w_k)
randomWalk = cumsum(DATA, 1);

% Gauss-Markov process of 1st order x_(k+1) = exp(-1/T_)
corrTimes = [2000, 500];

gaussMarkov = zeros([size(DATA),length(corrTimes)]);
gaussMarkov(1,:,1) = DATA(1,:); gaussMarkov(1,:,2) = DATA(1,:);

for kk = 1:length(corrTimes)
    for jj = 2:size(DATA,1)
        gaussMarkov(jj,:,kk) = exp(-1/corrTimes(kk))*gaussMarkov(jj-1,:,kk) + DATA(jj-1,:);
    end

    figure('Position',[0,0,800,1200]); subplot(4,1,1)  % Plot results
    set(groot,'DefaultAxesFontSize',14)
    set(groot,'DefaultLineLineWidth',1.5)
    for i=1:3
        plot(gaussMarkov(:,i,kk),cols(i)); hold on;
    end
    xlabel('Sample Number []'); ylabel('Signal value []')


    subplot(4,1,2)
    %fonts: set(groot,'DefaultAxesFontSize',14)
    %set(groot,'DefaultLineLineWidth',1.5)
    yRange = [0,0];
    for i = 1:3
        % Autocorrelation Using User defined Formulla
        phi = xcorr(gaussMarkov(:,i),'unbiased');
        %phi = autocorrelation(gaussMarkov(:,i,kk));
        Nc = length(phi);
        plot((1:Nc)-1,phi, cols(i)); hold on;
        yRange(1) = min(yRange(1), min(phi));
        yRange(2) = max(yRange(2), max(phi));
    end
    ylim([yRange(1)-10, yRange(2)+10]);
    xlabel('Sample Number []'); ylabel('Autocorrelation []')

    subplot(4,1,3)
    %fonts: set(groot,'DefaultAxesFontSize',14)
    %set(groot,'DefaultLineLineWidth',1.5)
    yRange = [1,1];

    Fs = 1; % Default frequency: 1Hz
    for i = 1:3 
        % Autocorrelation Using User defined Formulla
        [psd_randWalk1, freq1] = powerSpectralDensity(gaussMarkov(:,i,kk), Fs);
        loglog(freq1, psd_randWalk1, cols(i)); hold on;
        yRange(1) = min(yRange(1), min(psd_randWalk1));
        yRange(2) = max(yRange(2), max(psd_randWalk1));
    end
    xlabel('Sample Number [log]'); ylabel('Power Spectral Density [log]')
    xlim([min(freq1),max(freq1)]); ylim([yRange(1)*0.9, yRange(2)*1.1]);

    print(strcat('fig/gaussMarkov',num2str(corrTimes(kk))),'-dpdf')
    print(strcat('fig/gaussMarkov',num2str(corrTimes(kk))),'-depsc')

end


% %%
% 
% 
% % Plot Random Walk 
% figure('Position', [0 0 1200 600])
% plot(randomWalk(:,1),'b'); 
% hold on; %grid on;
% plot(randomWalk(:,2),'r')
% plot(randomWalk(:,3),'g')
% title('Random walk')
% xlabel('Sample number []'); ylabel('Signal value []')
% 
% % Gauss-Markov process of 1st order x_(k+1) = exp(-1/T_)
% corrTimes = [2000, 500];
% 
% gaussMarkov = zeros([size(whiteNoise),length(corrTimes)]);
% gaussMarkov(1,:,1) = whiteNoise(1,:); gaussMarkov(1,:,2) = whiteNoise(1,:);
% 
% for kk = 1:length(corrTimes)
%     for jj = 2:size(whiteNoise,1)
%         gaussMarkov(jj,:,kk) = exp(-1/corrTimes(kk))*gaussMarkov(jj-1,:,kk) + whiteNoise(jj-1,:);
%     end
% end
% 
% % Plot Random Walk 
% figure('Position', [0 0 1200 600])
% plot(gaussMarkov(:,1),'b'); 
% hold on; %grid on;
% plot(gaussMarkov(:,2),'r')
% plot(gaussMarkov(:,3),'g')
% title('gGauss-Markov')
% xlabel('Sample number []'); ylabel('Signal value []')
% 
% 
% % Plot in one plot
% figure('Position', [0 0 1200 600])
% plot(whiteNoise(:,1),'Color',[0,0,1]); 
% hold on; %grid on;
% plot(randomWalk(:,1),'Color',[0,1,0])
% plot(gaussMarkov(:,1,1),'Color',[0.5,0.5,0])
% plot(gaussMarkov(:,1,2),'Color',[0.6,0.5,1])
% title('Process for one random sequence')
% legend('White noise', 'Random walk', 'Gauss-Markov (T=2000)','Gauss-Markov (T=500)');
% xlabel('[]'); ylabel('Signal value []')
% 
% %% Exercise 2
% % Tableau noir
% % ac = xcorr(x, '');
% % plot(ac)
% % pwelch;
% % x_(k+1) = exp(-beta*Delta_t)*x_k+w_k
% % w_k~N(0,q)
% % q= sigma_GM ^2 *(1-exp(-2*beta*Detla_t));
% % x_(k+1) = phi*x_k + w_k
% % beta = 1/T
% % Delta_t = 1
% 
% close all;
% 
% % ------------------------------------
% %           Autocorleation
% % ------------------------------------
% 
% % Autocorrelation Using User defined Formulla
% phi = autocorrelation(whiteNoise(:,1));
% Nc = length(phi);
% figure; subplot(3,1,1) % Plot results
% plot((1:Nc)-1,phi,'b'); hold on;
% ylim([-0.2,1.2])
% 
% % Autocorrelation Using MATLAB defined Formula
% phi_x = xcorr(whiteNoise(:,1),'unbiased'); % unbiased autocorrelation estimate
% subplot(3,1,2)
% Nx = length(phi_x);
% phi_x = phi_x((Nx+1)/2:end);
% plot(1:(Nx+1)/2,phi_x,'b'); % Plot results
% ylim([-0.2,1.2])
% 
% % Plot error betwen the two functions
% err_corr = [phi;0]-phi_x; % add element in xcorr, as it's too short
% subplot(3,1,3); 
% plot(1:(Nx+1)/2,err_corr,'r');
% ylim([-1.2,1.2])
% 
% % Distribution meam square of error
% sum_err_low = mean(err_corr(1:1900).^2) 
% sum_err_high = mean(err_corr(1901:end).^2)
% 
% 
% %% ------------------------------------
% %         Power Spectral Density
% % ------------------------------------
% 
% Fs = 1; % Frequence: 1Hz
% [psd_whiteNoise1, freq1] = powerSpectralDensity(whiteNoise(:,1), Fs);
% figure; subplot(3,1,1)
% loglog(freq1, psd_whiteNoise1)
% xlim([min(freq1),max(freq1)])
% 
% [psd_randomWalk1, freq1] = powerSpectralDensity(randomWalk(:,1), Fs);
% subplot(3,1,2);
% loglog(freq1, psd_randomWalk1)
% xlim([min(freq1),max(freq1)])
% 
% [psd_gaussMarkov1, freq1] = powerSpectralDensity(gaussMarkov(:,1), Fs);
% subplot(3,1,3);
% loglog(freq1, psd_gaussMarkov1)
% xlim([min(freq1),max(freq1)])
% 
% 
% %% ------------------------------------
% %         Power Spectral Density
% % ------------------------------------

%%
fprintf('Program fnished. \n');