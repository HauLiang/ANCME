%% ---------------- Non-crossed modes Demo for ANCME -----------------------
%
% This is a simple example to test the ANCME algorithm 
% -- Adaptive Nonlinear Chirp Mode Estimation
%
% Author: Hao Liang 
%
% Last modified by: 22/12/09
%

clc; clear; close all

T = 1;     % time duration
SampFreq = 2000;  % sample frequency
t = 0:1/SampFreq:T;  % time variables

% Instantaneous amplitudes (IAs) 
a1 = sawtooth(9*pi*t,0.1)+2;
a2 = sawtooth(9*pi*t,0.1)+4;
a3 = sawtooth(9*pi*t,0.1)+6;

% Instantaneous frequencies (IFs)
f1 = 100 + 300*t + 9*pi*cos(9*pi*t);
f2 = 300 + 300*t + 9*pi*cos(9*pi*t);
f3 = 500 + 300*t + 9*pi*cos(9*pi*t);

% A three-component simulated nonlinear chirp signal (NCS)
y1 = a1.*cos(2*pi*(100*t+150*t.^2+sin(9*pi*t)));
y2 = a2.*cos(2*pi*(300*t+150*t.^2+sin(9*pi*t)));
y3 = a3.*cos(2*pi*(500*t+150*t.^2+sin(9*pi*t)));
y = y1 + y2 + y3;

% Show the signal's time-domain waveform and IFs
figure
set(gcf,'Position',[20 100 640 500]);
set(gcf,'Color','w');
plot(t,y,'b','linewidth',2);
xlabel('Time (s)','FontSize',24,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'linewidth',2);

figure
set(gcf,'Position',[20 100 640 500]);
set(gcf,'Color','w');
plot(t,[f1;f2;f3],'b','linewidth',2);
xlabel('Time (s)','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'linewidth',2);
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
ylim([80,820])


%% ANCME
beta = 1e-7;     % filter parameter
tol = 1e-5;      % tolerance of convergence criterion
iniIF = [100 + 300*t;300 + 300*t;500 + 300*t];   % initial IFs

% Start ANCME algorithm
tic;
[estIF, estIA, estMode] = ANCME(y, SampFreq, iniIF, beta, tol);
toc


%% The SNRs and REs of the estimated modes and frequencies 
SNR1 =  20*log10(norm(y1)/norm(y1 - estMode(1,:,end)))
SNR2 =  20*log10(norm(y2)/norm(y2 - estMode(2,:,end)))
SNR3 =  20*log10(norm(y3)/norm(y3 - estMode(3,:,end)))

RE1 =  norm(estIF(1,:,end)-f1)/norm(f1)
RE2 =  norm(estIF(2,:,end)-f2)/norm(f2)
RE3 =  norm(estIF(3,:,end)-f3)/norm(f3)


%% Show the estimated IF
figure
plot(t,[f1;f2;f3],'b','linewidth',2); hold on;
plot(t,estIF(:,:,end),'r','linewidth',2); 
set(gcf,'Position',[20 100 640 500]);	 
xlabel('Time (s)','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',24,'FontName','Times New Roman');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);
set(gcf,'Color','w');	
ylim([50,850])

