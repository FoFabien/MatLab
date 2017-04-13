function sd=modulateur_v2(M, phase_ini)

%Variable
Fs = 1;
Nb = 100;
N = Nb / log2(M);
Db = 0.1;
D = Db / log2(M);
T = 1 / D;
NPTS = 1000;
NFFT = 2 ^nextpow2(NPTS);
f = linspace(-Fs/2,Fs/2,NFFT);
fc = 0.2;

% Calcul
alpha_k=randint(ceil(Nb / nextpow2(M)),nextpow2(M)); % on organise les bits sur un nombre de colonne égale à la puissance de 2 de M. La fonction reshape pourrait aussi faire l'affaire
entiers = bi2de(alpha_k); % afin qu'ils soient convertis par bi2de

symbole=pskmod(entiers, M, phase_ini); % puis on module
scatterplot(symbole);
sd = symbole;

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

% Affichage
figure;
subplot(3,1,1);
plot(real(s),'-k');
hold on;
subplot(3,1,2);
plot(real(a_filtre),'-');
subplot(3,1,3);
plot(imag(a_filtre),'-r');

figure;
% Densité Spectrale
% peigne Dirac
dspdirac = abs(fft(a, 2^nextpow2(NPTS))).^2 /NFFT;
subplot(3,1,1);
semilogy(f, fftshift(dspdirac));

% passe-bas
dsplowpass = abs(fft(a_filtre, 2^nextpow2(NPTS))).^2 /NFFT;
subplot(3,1,2);
semilogy(f, fftshift(dsplowpass));

% passe-bande
dspbandpass = abs(fft(s, 2^nextpow2(NPTS))).^2 /NFFT;
subplot(3,1,3);
semilogy(f, fftshift(dspbandpass));