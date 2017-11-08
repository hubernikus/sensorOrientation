%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               Sensor Orientation - EPFL
%               Lab 2 - 
%               Lukas Huber - 2017/10/13
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accelerometer in x & z direction
% Sensors used: 3 - FSAS/IMAR (tactical, FOG) 6 - Smart-phone (low-cost,
% MEMS) 
clear all;

%%

clc; clear variables; close all;

addpath(genpath('../Ex1'))
addpath(genpath('../functions'))

dataPath = 'imuData2017/'; 


%% Define General Parameters 

c = ['r','g','b']; 


%% Notes Black-Board
%
% d = [t, x_0 y_0 z_0, x_A y_A z_A]
% x_k = exp(-beta Dt) * x_{k-1} + w_k (NOT \simga_w ?)
% (simga_GM)^2 = (exp(-beta Dt)^2 * \sigma_GM - (sigma_W)^2
% sigma_W^2 = sigma _{GM1}^2 - exp(-2 beta Dt) sigma_GM1^2
% simga_w = sqrt(1-exp(-2 beta Dt) sigma_GM^2
% w_k = sqrt(sigma_W^2) * randn 


%% Import Data
[dataIMAR, fIMU] = readimu(strcat(dataPath, 'imu3_20171006_imar-fsas.imu'),'IMAR');
load(strcat(dataPath, 'imu6_Phone.mat'));

%% IMAR Analysis
%N_sample = 'end';
N_sampleIMAR = size(dataIMAR,1);


dt_IMAR = round(mean(dataIMAR(2:N_sampleIMAR,1) - dataIMAR(1:N_sampleIMAR-1,1)),5);

figure;
for ii = [1,3]
    plot(dataIMAR(:,1)-dataIMAR(1,1), dataIMAR(:,ii+4),c(ii)); hold on;
end
xlim([dataIMAR(1,1),dataIMAR(N_sampleIMAR,1)]-dataIMAR(1,1))
grid on;
xlabel('Time [s]','Interpreter','latex'); ylabel('Accelerometer measurement [$m/s^2$]','Interpreter','latex')
legend('x direction','z direction')
print(strcat('fig/','IMAR_xz'),'-depsc')
%% Phone Analysis
N_samplePhone = size(t_a,1);

dt_phone = round(mean(t_a(2:N_samplePhone)-t_a(1:N_samplePhone-1)),5);

figure(2);
for ii = [1,3]
    plot(t_a(1:N_samplePhone), a(:,ii),c(ii)); hold on;
end
xlim([t_a(1,1),t_a(N_samplePhone,1)])

xlabel('Time [s]','Interpreter','latex'); 
ylabel('Accelerometer measurement [$m/s^2$]','Interpreter','latex')
legend('x direction','z direction')
grid on;
print(strcat('fig/','phone_xz'),'-depsc')

%% Correlation analysis
%figure(3);
DATA = [dataIMAR(:,6),dataIMAR(:,4)];
[ sigma_corr_imar, sigma_sig_imar, ind_e1_imar ] = signalAnalysis(DATA,'Datanalysis_IMAR',dt_IMAR)

T_imar = ind_e1_imar*dt_IMAR
beta_imar = 1/ind_e1_imar*dt_IMAR; 


%% z first, x second 
DATA = [a(:,3),a(:,1)];
[ sigma_corr_phone, sigma_sig_phone, ind_e1_phone ] = signalAnalysis(DATA,'Datanalysis_mobilePhone',dt_phone)

T_phone = ind_e1_phone*dt_phone
beta_phone = 1/ind_e1_phone*dt_phone;

%% Allan deviation
allanDev = 0; % 1 to plot the allan deviation, very timeconsuming...
if(allanDev)
    allan_DATA = allandev(a(:,3),'AllanDev_phoneZ');
    print(strcat('fig/','allanDev_phoneZ'),'-depsc')

    allan_DATA = allandev(a(:,1),'AllanDev_phoneX');
    print(strcat('fig/','allanDev_phoneX'),'-depsc')

    allan_DATA = allandev(dataIMAR(:,4),'AllanDev_dataImarZ');
    print(strcat('fig/','allanDev_dataImarZ'),'-depsc')

    allan_DATA = allandev(dataIMAR(:,2),'AllanDev_dataImarX');
    print(strcat('fig/','allanDev_dataImarX'),'-depsc')
end


%% Print to file
printToFile = 0;
if(printToFile)
    fileID = fopen('phoneZ.txt','w');
    fprintf(fileID,'%4.6f \n',a(:,1)-mean(a(:,1)));
    fclose(fileID);

    fileID = fopen('phoneX.txt','w');
    fprintf(fileID,'%4.6f \n',a(:,3)-mean(a(:,3)));
    fclose(fileID);

    fileID = fopen('imarZ.txt','w');
    fprintf(fileID,'%4.6f \n',dataIMAR(:,4)-mean(dataIMAR(:,4)));
    fclose(fileID);

    fileID = fopen('imarX.txt','w');
    fprintf(fileID,'%4.6f \n',dataIMAR(:,2)-mean(dataIMAR(:,2)));
    fclose(fileID);
end


%% Simulation phone
a_sim = whiteNoiseGen(N_samplePhone,1,1.16e-2);
a_sim(:,2) = whiteNoiseGen(N_samplePhone,1,1.29e-2);

DATA = [a_sim(:,2),a_sim(:,1)];
[ sigma_corr_phone_sim, sigma_sig_phone_sim, ind_e1_phone_sim ] = signalAnalysis(DATA,'Datanalysis_mobilePhone_sim',dt_phone)

T_phone_sim = ind_e1_phone_sim*dt_phone
beta_phone_sim = 1/ind_e1_phone_sim*dt_phone;

%% Simulation IMAR
dataIMAR_simX = whiteNoiseGen(N_sampleIMAR,1,6.00e-4);
[ sigma_corr_imar_simx, sigma_sig_imar_simx, ind_e1_imar_simx ] = signalAnalysis(dataIMAR_simX,'Datanalysis_mobileImar_simX',dt_phone)


T_phone_simx = ind_e1_imar_simx*dt_phone
beta_phone_simx = 1/ind_e1_imar_simx*dt_phone;

%%
freqSin = 100*1.14;   % why multiply by this factor
magSin = sqrt((0.015)/(pi/2)*0.002);
dataIMAR_simZ = whiteNoiseGen(N_sampleIMAR,1,6.90e-3)+sin((1:N_sampleIMAR)'*freqSin)*magSin;
[sigma_corr_phone_simz, sigma_sig_phone_simz, ind_e1_imar_simz ] = signalAnalysis(dataIMAR_simZ,'Datanalysis_mobilePhone_simZ',dt_IMAR)

T_phone_simz = ind_e1_imar_simz*dt_phone
beta_phone_simz = 1/ind_e1_imar_simz*dt_phone;

%% Allan deviation
allanDev = 1; % 1 to plot the allan deviation, very timeconsuming...
if(allanDev)        
    allan_DATA = allandev(a_sim(:,2),'AllanDev_phoneZ_sim');
    print(strcat('fig/','allanDev_phoneZ_sim'),'-depsc')

    allan_DATA = allandev(a_sim(:,1),'AllanDev_phoneX_sim');
    print(strcat('fig/','allanDev_phoneX_sim'),'-depsc')

    allan_DATA = allandev(dataIMAR_simZ,'AllanDev_dataImarZ_sim');
    print(strcat('fig/','allanDev_dataImarZ_sim'),'-depsc')

    allan_DATA = allandev(dataIMAR_simX,'AllanDev_dataImarX_sim');
    print(strcat('fig/','allanDev_dataImarX_sim'),'-depsc')
end

