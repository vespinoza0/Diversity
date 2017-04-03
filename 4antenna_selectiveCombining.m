clear
close all
 
SNR_dB_range = [0:1:25];
 
M = 1;           % # of antennas
%M = 2;          % # of antennas
%M = 3;          % # of antennas
N = 2e6;         % # of samples
 
s1 = (randi([0 1],N,1)*2-1 + 1j*(randi([0 1],N,1)*2-1));   %QPSK1*sqrt(2) - using Eb/N0

h = (randn(N,M) + 1j*randn(N,M))/sqrt(2);    %Rayleigh channel_1
h1 = (randn(N,M) + 1j*randn(N,M))/sqrt(2);   %Rayleigh channel_2
h2 = (randn(N,M) + 1j*randn(N,M))/sqrt(2);   %Rayleigh channel_3
h3 = (randn(N,M) + 1j*randn(N,M))/sqrt(2);   %Rayleigh channel_4
noise = (randn(N,M) + 1j*randn(N,M))/sqrt(2);  
 
BER(length(SNR_dB_range),M)=0;  

s_rec1 = abs(h) .* s1;        %received signal; after phase synchronized
s_rec2 = abs(h1) .* s1;       %received signal; after phase synchronized
s_rec3 = abs(h2) .* s1;       %received signal; after phase synchronized
s_rec4 = abs(h3) .* s1;       %received signal; after phase synchronized
y1 = s_rec1;        %% 1 antenna
y2 = max(y1, s_rec2);  %%best received signal with 2 antennas
y3 = max(y2, s_rec3);  %%best received signal with 3 antennas
y4 = max(y3, s_rec4);  %%best received signal with 4 antennas
        
for iSNR = 1:length(SNR_dB_range)  %for i = 1:26
sigma = 10^(-SNR_dB_range(iSNR)/20); %%
z1 = y1 + sigma*noise(:,1);  % z =received signal with noise
ydet = sign(real(z1)) + 1j*sign(imag(z1)); %%ydet= decision =(+/-)1 +/- j; 
BER1(iSNR,M) = sum(abs(real(ydet - s1)) + abs(imag(ydet - s1)) )/N/2/2;
end

for iSNR = 1:length(SNR_dB_range)  %for i = 1:26 
sigma = 10^(-SNR_dB_range(iSNR)/20); %%
z2 = y2 + sigma*noise(:,1);  % y =received signal with noise
ydet = sign(real(z2)) + 1j*sign(imag(z2)); %%ydet= decision =(+/-)1 +/- j; 
BER2(iSNR,M) = sum(abs(real(ydet - s1)) + abs(imag(ydet - s1)) )/N/2/2;
end

for iSNR = 1:length(SNR_dB_range)  %for i = 1:26 
sigma = 10^(-SNR_dB_range(iSNR)/20); %%
z3 = y3 + sigma*noise(:,1);  % y =received signal with noise
ydet = sign(real(z3)) + 1j*sign(imag(z3)); %%ydet= decision =(+/-)1 +/- j; 
BER3(iSNR,M) = sum(abs(real(ydet - s1)) + abs(imag(ydet - s1)) )/N/2/2;
end

for iSNR = 1:length(SNR_dB_range)  %for i = 1:26 
sigma = 10^(-SNR_dB_range(iSNR)/20); %%
z4 = y4 + sigma*noise(:,1);  % y =received signal with noise
ydet = sign(real(z4)) + 1j*sign(imag(z4)); %%ydet= decision =(+/-)1 +/- j; 
BER4(iSNR,M) = sum(abs(real(ydet - s1)) + abs(imag(ydet - s1)) )/N/2/2;
end

a1 = BER1(:,1);
a2 = BER2(:,1);
a3 = BER3(:,1);
a4 = BER4(:,1);
b = SNR_dB_range;

figure(2)
semilogy(b,a1,b,a2,b,a3,b,a4);
grid on
legend('M=1','M=2','M=3','M=4');
xlabel('E_b/N_0 (dB)')
ylabel('BER')
axis([0 25 .9999e-5 1])
