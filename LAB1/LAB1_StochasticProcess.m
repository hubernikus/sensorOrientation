%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      EPFL - Sensor Orientation
%                     LAB 1 - Stochastic Processes
%                              Lukas Huber
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear variables; close all; 

% Initialize (pseudo-)random number generator based on birthday
rng(19920526);  

fprintf('Program started. \n');

%% Part A

N_random = 200000;  % Length of sequences

%% Generation of 3 random white noise sequences
whiteNoise = whiteNoiseGen(N_random);
[std_WN] = signalAnalysis(whiteNoise, 'whiteNoise')


fileID = fopen('whiteNoise.txt','w');
fprintf(fileID,'%4.6f \n',whiteNoise(:,1));
fclose(fileID);

%% Generation of 3 random random walk
randomWalk = randWalkGen(N_random);
[std_WN] = signalAnalysis(randomWalk, 'randomWalk')

fileID = fopen('randomWalk.txt','w');
fprintf(fileID,'%4.6f \n',randomWalk(:,1));
fclose(fileID);

%% Generation of 3 gauss markov sequences
corrTime = 2000;
gaussMarkov = gaussMarkovGen(N_random, corrTime);
[std_WN] = signalAnalysis(gaussMarkov, 'gaussMarkov2000')

fileID = fopen('whiteNoise2000.txt','w');
fprintf(fileID,'%4.6f \n',whiteNoise(:,1));
fclose(fileID);

%% Generation of 3 gauss markov sequences
corrTime = 500;
gaussMarkov = gaussMarkovGen(N_random, corrTime);
[std_WN] = signalAnalysis(gaussMarkov, 'gaussMarkov500')

fileID = fopen('gaussMarkov500.txt','w');
fprintf(fileID,'%4.6f \n',gaussMarkov(:,1));
fclose(fileID);

%%
fprintf('Program fimished. \n');