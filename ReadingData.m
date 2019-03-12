%Reading Data from ICA algorithm

clear 
clc

tic

%Retrieving data from textfile
ICAsourceData = fopen('SourcesEstimation.txt');
N=fscanf(ICAsourceData,'%d',[1 1]);
M=fscanf(ICAsourceData,'%g',[1 1]);
Sest=fscanf(ICAsourceData,'%g %g',[N M]);
timeVector=fscanf(ICAsourceData,'%g %g',[1 M]);

fclose(ICAsourceData);

figure
plot(timeVector, Sest(1,:))
xlabel('time (s)') 
ylabel('Signal Amplitude') 
legend('Source Estimation 1')

figure
plot(timeVector, Sest(2,:))
xlabel('time (s)') 
ylabel('Signal Amplitude') 
legend('Source Estimation 2')

figure
plot(timeVector, Sest(3,:))
xlabel('time (s)')
ylabel('Signal Amplitude') 
legend('Source Estimation 3')

figure
plot(timeVector, Sest(4,:))
xlabel('time (s)')
ylabel('Signal Amplitude') 
legend('Source Estimation 4')

figure
plot(timeVector, Sest(5,:))
xlabel('time (s)')
ylabel('Signal Amplitude') 
legend('Source Estimation 5')

figure
plot(timeVector, Sest(6,:))
xlabel('time (s)')
ylabel('Signal Amplitude') 
legend('Source Estimation 6')

toc