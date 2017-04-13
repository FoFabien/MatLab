%function sd=demodulateur(M, signal_recu)
clc;
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
fc2 = fc * 0.01;
delta = 1 / (2*T);
phi = 0;
phi2 = pi / 8;

% création d'un signal BPSK. On peut aussi utilise modulateur.m
alpha_k=randint(ceil(Nb / nextpow2(M)),nextpow2(M));
entiers = bi2de(alpha_k);
s_mod=pskmod(entiers, M);
porte=ones(T,1);
s_lowpass=filter(porte,T,s_mod); % signal passe-bas
for(t1=0:N-1)
signal_recu(t1+1,1) = s_lowpass(t1+1,1) .* cos(2 * pi * fc * t1); % signal passe-bande
signal_recu2(t1+1,1) = s_lowpass(t1+1,1) .* cos(2 * pi * fc2 * t1); % signal passe-bande AVEC fc+1%
end;


% création des deux porteuses pures et retour en bande de base sur 2 voies
for(t1=0:N-1)
s1(t1+1,1)= signal_recu(t1+1,1) * cos(2 * pi * fc * t1 + phi);
s2(t1+1,1)= - signal_recu(t1+1,1) * sin(2 * pi * fc * t1 + phi);
end;

for(t1=0:N-1) % même chose mais avec fc+1%
s3(t1+1,1)= signal_recu(t1+1,1) * cos(2 * pi * fc2 * t1 + phi2);
s4(t1+1,1)= - signal_recu(t1+1,1) * sin(2 * pi * fc2 * t1 + phi2);
end;

% filtrage adapté
adapt_porte = fliplr(porte); % création du filtrage adapté au filtre de mise en forme
adapt1 = filter(adapt_porte, 1, s1); % et on filtre les deux voies
adapt2 = filter(adapt_porte, 1, s2);
% même chose pour fc+1%
adapt3 = filter(adapt_porte, 1, s3); 
adapt4 = filter(adapt_porte, 1, s4);

% échantillonne au rythme symbole T
s_echantillon1 = adapt1(T:T:N,1);
s_echantillon2 = adapt2(T:T:N,1);

s_final = adapt1 + 1i * adapt2; % somme des deux voies
s_final2 = adapt3 + 1i * adapt4; % même chose pour fc+1%

s_demod = pskdemod(s_final, M); % modulation
s_demod2 = pskdemod(s_final2, M); % même chose pour fc+1%

%représentation temporelle -------------------------------------------
tt = 0.1:0.1:10;

s = signal_recu(T:T:N,1);

figure
subplot(2,1,1);
plot(tt(T:T:N), real(s), 'ok');
hold on;
plot(tt, real(s_lowpass));
hold on;
plot(tt,real(adapt1 + 1i * adapt2), '-r');
hold on;
plot(tt(T:T:N), real(s_echantillon1 + 1i * s_echantillon2), 'og');

subplot(2,1,2);
plot(tt(T:T:N), imag(s), 'ok');
hold on;
plot(tt, imag(s_lowpass));
hold on;
plot(tt,imag(adapt1 + 1i * adapt2), '-r');
hold on;
plot(tt(T:T:N), imag(s_echantillon1 + 1i * s_echantillon2), 'og');

%représentation fréquentielle -------------------------------------------
figure
subplot(2,1,1)
semilogy(f, fftshift(abs(fft(real(adapt1 + 1i * adapt2), 2^nextpow2(npts))).^2 /NFFT));
hold on
semilogy(f, fftshift(abs(fft(real(adapt3 + 1i * adapt4), 2^nextpow2(npts))).^2 /NFFT), '-r'); % affichage pour fc+1%
hold on;
subplot(2,1,2)
semilogy(f, fftshift(abs(fft(signal_recu, 2^nextpow2(npts))).^2 /NFFT));
hold on;
semilogy(f, fftshift(abs(fft(signal_recu2, 2^nextpow2(npts))).^2 /NFFT), '-r'); % affichage pour fc+1%