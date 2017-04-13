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



sum_symbole_fsk = zeros(NFFT,1); % même chose que pour modulateur_v1.m, ces signaux seront utilisés pour la moyenne des estimations
sum_signal_fsk = zeros(1, NFFT);

for (jj=1:Nstat);
% Calcul
alpha_k=randint(ceil(Nb / nextpow2(M)),nextpow2(M));
entiers = bi2de(alpha_k);

% modulation avec FSKMOD
symbole_fsk=fskmod(entiers, M, delta, T, Fs,'discont'); % on module en discontinuité de phase avec fskmod

e = zeros(N*T, 1);
for(ii=1:N*T)
   e(ii,1) =  exp(j * 2 * pi * fc * (ii - 1)); % on créé un signal exponentiel
end;

t1=0:npts-1;
% modulation sans FSKMOD (discontinuité de phase)
c1 = cos(2 * pi * (fc - delta / 2) * t1); % on prépare les deux cosinus utilisés
c2 = cos(2 * pi * (fc + delta / 2) * t1);

%signal_fsk = zeros(N*T,1);
for(ii=1:Nb)
        if(alpha_k(ii) == 0) % si le bit est nul
            signal_fsk((ii-1)*T+1:ii*T) = c1((ii-1)*T+1:ii*T); % on prend la valeur de cos(2pi (fc-delta/2))
        else
            signal_fsk((ii-1)*T+1:ii*T) = c2((ii-1)*T+1:ii*T); % sinon on prend la valeur de cos(2pi (fc+delta/2))
        end
end;

sum_symbole_fsk = sum_symbole_fsk + abs(fft(real(symbole_fsk .* e), 2^nextpow2(npts))).^2 /NFFT; % somme des fft comme dans modulateur_v1.m
sum_signal_fsk = sum_signal_fsk + abs(fft(signal_fsk, 2^nextpow2(npts))).^2 /NFFT;

end;

%représentation temporelle
figure
subplot(2, 1, 1);
plot(real(symbole_fsk .* e)); % on multiplie par un exponentiel pour avoir le signal réel
subplot(2, 1, 2);
plot(reshape(signal_fsk,1,Nb*T), '-r'); % on affiche ici le signal obtenue SANS fskmod

%représentation fréquentielle
figure
semilogy(f, fftshift(sum_symbole_fsk) / Nstat, '-rs'); % et ici on affiche les représentations fréquencielles
hold on;
semilogy(f, fftshift(sum_signal_fsk) / Nstat);
hold on;