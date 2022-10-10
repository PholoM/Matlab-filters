% task 1: to read and understand the audio file in matlab
% task 2: to develop lowpass
% task 3: highpass
% task 4: bandpass
% task 5: bandstop

clc, clear, close;

fc = 3000; %cut off frequency in Hz
[s1, fs] = audioread('white_noise.wav');
N = length(s1);
Ts = 1/fs;
t = 0:Ts:(N-1)*Ts;
Fc  = fc/fs;        %normalised cut off frequency

transition = 1000;
TRx = transition/fs;
dF_rec = 0.9;   % 0.9/N
% dF_han = 3.1;
% dF_ham = 3.3;
% dF_bla = 5.5;


%low pass FIR filter
Nf_rec = dF_rec/ TRx; 
% Nf_han = dF_han/ Fc;
% Nf_ham = dF_ham/ Fc; 
% Nf_bla = dF_bla/ Fc; 

Nf = floor(Nf_rec);
if (mod(Nf, 2) ==0)
    Nf = Nf+1;
end

% M_rec =  (Nf_rec - 1)/2; % number of truncation samples
% M_han =  (Nf_han - 1)/2;
M = (Nf - 1)/2;


%Nf = 2*M+1;      %total length of impulse response truncated
hf = zeros(1, M);

for n = 1:M
    hf(n) = 2*Fc*sin(2*pi*Fc*n)/(n*2*pi*Fc); %irt
end
%hf = 2*Fc*sin(n*2*pi*Fc)/(n*2*pi*Fc);   


hb = fliplr(hf);
hd = [hb 2*Fc hf];

% incorporate window here
% 1. hanning window whn
whn = zeros(1, Nf);
for n = 1:Nf
    whn(n) = 0.5 - cos((2*pi*n)/(N-1));
end
hwhn = hd.*whn;

% 2. hamming window whm
whm = zeros(1, Nf);
for n = 1:Nf
    whm(n) = 0.54 - 0.46*cos((2*pi*n)/(N-1));
end
hwhm = hd.*whm;

% 3. Blackmann window wbl
wbl = zeros(1, Nf);
for n = 1:Nf
    wbl(n) = 0.42 - 0.5*cos((2*pi*n)/(N-1)) +  0.08*cos((2*pi*n)/(N-1));
end
hwbl = hd.*wbl;

figure;
plot(1:Nf, hd);
hold on;
plot(1:Nf, hwhn);
plot(1:Nf, hwhm);
plot(1:Nf, hwbl);
title("Low-Pass Impulse Response with windows");
legend("Rectangular", "Hanning", "Hamming", "Blackman");
hold off;
pause;




s1f = conv(s1, hd);
s1hn = conv(s1, hwhn);
s1hm = conv(s1, hwhm);
s1bl = conv(s1, hwbl);

%sound(s1, fs);
pspectrum(s1, fs);
%pause;
hold on;

pspectrum(s1f, fs);
%sound(s1f, fs);
%pause;


pspectrum(s1hn, fs);
%sound(s1hn, fs);
%pause;

pspectrum(s1hm, fs);
%sound(s1hm, fs);
%pause;

pspectrum(s1bl, fs);
%sound(s1bl, fs);
%pause;
legend('unfiltered','rectangular', 'hanning', 'hamming', 'blackman');
%experiment

hold off


