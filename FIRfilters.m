% task 1: to read and understand the audio file in matlab
% task 2: to develop lowpass
% task 3: highpass
% task 4: bandpass
% task 5: bandstop

clc, clear, close;

% fc = 50; %cut off frequency in Hz
% [s1, fs] = audioread('audio.wav');
% N = length(s1);
% Ts = 1/fs;
% t = 0:Ts:(N-1)*Ts;
% Fc  = fc/fs;        %normalised cut off frequency


%low pass FIR filter
fs = 44100;
fc = 2500; %cut off frequency in Hz
Fc  = fc/fs;        %normalised cut off frequency
M = 341; % number of truncation samples
Nf = 2*M+1;      %total length of impulse response truncated
hf = zeros(1, M);

for n = 1:M
    hf(n) = 2*Fc*sin(2*pi*Fc*n)/(n*2*pi*Fc); %irt
end
%hf = 2*Fc*sin(n*2*pi*Fc)/(n*2*pi*Fc);   
hb = fliplr(hf);
hd = [hb 2*Fc hf];



