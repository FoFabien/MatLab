clc

close all;
clear all;

%Variable
Fs = 1;
M = 2;
Nb = 100; % 10 en temporel
N = Nb / log2(M);
Db = 0.1;
D = Db / log2(M);
T = 1 / D;
npts = N*T;
NFFT = 2 ^nextpow2(npts);
f = linspace(-Fs/2,Fs/2,NFFT);
fc = 0.2;
delta = 1 / (2*T);
Nstat = 100;

% même chose que dans modulateur_v4_disc.m

sum_fsk_discont = zeros(NFFT,1);
sum_fsk_cont = zeros(NFFT,1);

for (jj=1:Nstat);
% Calcul
alpha_k=randint(ceil(Nb / nextpow2(M)),nextpow2(M));
entiers = bi2de(alpha_k);

% modulation avec FSKMOD
fsk_discont=fskmod(entiers, M, delta, T, Fs,'discont'); % on module avec FSKMOD en discontinuité de phase
fsk_cont=fskmod(entiers, M, delta, T, Fs); % et en continuité de phase

e = zeros(N*T, 1);
for(ii=1:N*T)
   e(ii,1) =  exp(j * 2 * pi * fc * (ii - 1));
end;

sum_fsk_discont = sum_fsk_discont + abs(fft(real(fsk_discont .* e), 2^nextpow2(npts))).^2 /NFFT;
sum_fsk_cont = sum_fsk_cont + abs(fft(real(fsk_cont .* e), 2^nextpow2(npts))).^2 /NFFT;

end;

%représentation temporelle
figure
subplot(2, 1, 1);
plot(real(fsk_discont .* e));
subplot(2, 1, 2);
plot(real(fsk_cont .* e), '-r');

%représentation fréquentielle
figure
semilogy(f, fftshift(sum_fsk_discont) / Nstat);
hold on;
semilogy(f, fftshift(sum_fsk_cont) / Nstat, '-r');
hold on;