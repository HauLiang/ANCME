%% ---------------- Crossed-modes Demo for ANCME ---------------------------
%
% This is a simple example to test the ANCME algorithm 
% -- Adaptive Nonlinear Chirp Mode Estimation
%
% Author: Hao Liang 
%
% Last modified by: 22/12/09
%

clc; clear; close all

T = 15;     % time duration
SampFreq = 2000;  % sample frequency
t = 0:1/SampFreq:1;  % time variables

% Instantaneous amplitudes (IAs) and instantaneous frequencies (IFs)
a4 = exp(-0.06*t)+t+0.5*sawtooth(25*pi*t,1);
f4 = 100 + 300*t + 9*pi*cos(9*pi*t);
a5 = exp(-0.03*t)+t+0.5*sawtooth(20*pi*t,1);
f5 = 400 - 300*t + 9*pi*cos(9*pi*t);

% A two-component simulated nonlinear chirp signal (NCS)
y4 = a4.*cos(2*pi*(15+100*t+150*t.^2+sin(9*pi*t)));
y5 = a5.*cos(2*pi*(20+400*t-150*t.^2+sin(9*pi*t)));
g = y4 + y5;

% Show the signal's time-domain waveform and IFs
figure
set(gcf,'Position',[20 100 640 500]);
set(gcf,'Color','w');
plot(t,g,'b','linewidth',2);
xlabel('Time (s)','FontSize',24,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',24,'FontName','Times New Roman');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'FontSize',24);
set(gca,'linewidth',2);

figure
set(gcf,'Position',[20 100 640 500]);
set(gcf,'Color','w');
plot(t,[f4;f5],'b','linewidth',2);
xlabel('Time (s)','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',24,'FontName','Times New Roman');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'FontSize',24);
set(gca,'linewidth',2);


%% ANCME
beta = 1e-8;    % filter parameter
tol = 1e-3;     % tolerance of convergence criterion
iniIF = [100+300*t; 400-300*t];     % initial IFs

% Start ANCME algorithm
tic
[estIF, estIA, estMode] = ANCME(g, SampFreq, iniIF, beta, tol);
toc


%% The SNRs and REs of the estimated modes and frequencies 
SNR1 =  20*log10(norm(y4)/norm(y4 - estMode(1,:,end)))
SNR2 =  20*log10(norm(y5)/norm(y5 - estMode(2,:,end)))

RE1 =  norm(estIF(1,:,end)-f4)/norm(f4)
RE2 =  norm(estIF(2,:,end)-f5)/norm(f5)


%% Show the estimated IF
figure; x1 = 0.4; y4 = 180; x2 = 0.6; y5 = 320;
plot(t,[f4;f5],'b','linewidth',1.5);
hold on;
plot(t,estIF(1,:,end),'r','linewidth',1.5); hold on;
plot(t,estIF(2,:,end),'m','linewidth',1.5);
set(gcf,'Position',[20 100 640 500]);
set(gcf,'Color','w');
xlabel('Time (s)','FontName','Times New Roman');
ylabel('Frequency (Hz)','FontName','Times New Roman');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
ylim([0 500])
set(gca,'FontSize',24)
set(gca,'linewidth',2);
rectangle('Position',[x1 y4 x2-x1 y5-y4],'EdgeColor','k','Linewidth',1);
h1 = axes('position',[0.62 0.20 0.25 0.25]);
axis(h1);
plot(t,[f4;f5],'b','linewidth',1.5);
hold on;
plot(t,estIF(1,:,end),'r','linewidth',1.5); hold on;
plot(t,estIF(2,:,end),'m','linewidth',1.5);
xlim([x1 x2]);ylim([y4 y5]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
set(gca,'fontsize',12)

