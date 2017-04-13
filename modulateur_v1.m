close all;
clear all;

%Variable
Fs = 1;
M = 2;
Nb = 100;
N = Nb / log2(M);
Db = 0.1;
D = Db / log2(M);
T = 1 / D;
NPTS = 1000;
NFFT = 2 ^nextpow2(NPTS);
f = linspace(-Fs/2,Fs/2,NFFT);
fc = 0.2;
Nstat = 100; % nombre d'estimation


sumdspdirac = zeros(NFFT,1); % préparation des 3 signaux qui seront utilisés pour moyenné
sumdsplowpass = zeros(NFFT,1);
sumdspbandpass = zeros(NFFT,1);

for (ii=1:Nstat); % début de la boucle permettant de faire Nstat estimation

% Même chose que dans Modulateur_v0.m ----------------------------------------------
% Calcul
alpha_k=randint(Nb,1);
symbole=pskmod(alpha_k, M);

%Signal Symboles
a = zeros(N * T, 1);
for(i=0:N-1)
    a(i * T + 1) = symbole(i + 1);
end;

%Signal porte
porte=ones(T,1);

%Signal Passe-Bas
a_filtre=filter(porte,1,a);

%Signal Passe-Bande
c = zeros(N * T, 1);
for(i=1:N*T)
    c(i) = cos(2*pi*fc*i);
end;
s = a_filtre .* c;

% Densité Spectrale
% peigne Dirac
dspdirac = abs(fft(a, 2^nextpow2(NPTS))).^2 /NFFT;
subplot(3,1,1);
semilogy(f, fftshift(dspdirac));
hold on;

% passe-bas
dsplowpass = abs(fft(a_filtre, 2^nextpow2(NPTS))).^2 /NFFT;
subplot(3,1,2);
semilogy(f, fftshift(dsplowpass));
hold on;

% passe-bande
dspbandpass = abs(fft(s, 2^nextpow2(NPTS))).^2 /NFFT;
subplot(3,1,3);
semilogy(f, fftshift(dspbandpass));
hold on;
% ---------------------------------------------------------------------------------

sumdspdirac = dspdirac + sumdspdirac; % on ajoute les signaux
sumdsplowpass = dsplowpass + sumdsplowpass;
sumdspbandpass = dspbandpass + sumdspbandpass;

end;

subplot(3,1,1);
semilogy(f, fftshift(sumdspdirac) / Nstat, '-r'); % et on affiche la moyenne (ce qui reviens à diviser la somme des estimations par Nstat)
subplot(3,1,2);
semilogy(f, fftshift(sumdsplowpass) / Nstat, '-r');
subplot(3,1,3);
semilogy(f, fftshift(sumdspbandpass) / Nstat, '-r');