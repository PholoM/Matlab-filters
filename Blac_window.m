% task 1: to read and understand the audio file in matlab
% task 2: to develop lowpass
% task 3: highpass
% task 4: bandpass
% task 5: bandstop

clc, clear, close;

fc = 3000; %cut off frequency in Hz
f1 = 2000;
f2 = 4000;
[s1, fs] = audioread('white_noise.wav');
N = length(s1);
Ts = 1/fs;
t = 0:Ts:(N-1)*Ts;
Fc  = fc/fs;        %normalised cut off frequency
F1 = f1/fs;
F2 = f2/fs;

transition = 1000;
TRx = transition/fs;
dF_rec = 5.5;   % 0.9/N



%low pass FIR filter
Nf_rec = dF_rec/ TRx; 


Nf = floor(Nf_rec);
if (mod(Nf, 2) ==0)
    Nf = Nf+1;
end

M = (Nf - 1)/2;


%Nf = 2*M+1;      %total length of impulse response truncated
hf = zeros(1, M);

for n = 1:M
    hf(n) = 2*Fc*sin(2*pi*Fc*n)/(n*2*pi*Fc); %irt
end


hb = fliplr(hf);
hd = [hb 2*Fc hf];

wbl = zeros(1, Nf);
for n = 1:Nf
    wbl(n) = 0.42 - 0.5*cos((2*pi*n)/(N-1)) +  0.08*cos((2*pi*n)/(N-1));
end
hwbl = hd.*wbl;

lp_filtered = conv(s1, hwbl);

audiowrite("Blackman LP.wav", lp_filtered, 44100);

%% high pass
hpf = zeros(1, M);
for n = 1:M
    hpf(n) = 1 - (2*Fc*sin(2*pi*Fc*n)/(n*2*pi*Fc)); %irt
end

hpb = fliplr(hpf);
hpd = [hpb 1-2*Fc hpf];
hpwbl = hpd.*wbl;
hp_filtered = conv(s1, hpwbl);
audiowrite("Blackman HP.wav", hp_filtered, 44100);

%% band pass
bpf = zeros(1, M);
for n = 1:M
    bpf(n) = 2*F2*sin(2*pi*F2*n)/(n*2*pi*F2) - 2*F1*sin(2*pi*F1*n)/(n*2*pi*F1); %irt
end
%hf = 2*Fc*sin(n*2*pi*Fc)/(n*2*pi*Fc);   
bpb = fliplr(bpf);
bpd = [bpb 2*F2-2*F1 bpf];
bpwbl = bpd.*wbl;
bp_filtered = conv(s1, bpwbl);
audiowrite("Blackman BP.wav", bp_filtered, 44100);


%% band stop
bsf = zeros(1, M);
for n = 1:M
    bsf(n) = 1 - (2*F2*sin(2*pi*F2*n)/(n*2*pi*F2) - 2*F1*sin(2*pi*F1*n)/(n*2*pi*F1)); %irt
end

bsb = fliplr(bsf);
bsd = [bsb 1-(2*F2-2*F1) bsf];
bswbl = bsd.*wbl;
bs_filtered = conv(s1, bswbl);
audiowrite("Blackman BS.wav", bs_filtered, 44100);


