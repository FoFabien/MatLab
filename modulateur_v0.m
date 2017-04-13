close all;
clear all;

%Variable
Fs = 1; % fréquence d'échantillonnage
M = 2; 
Nb = 100; % nombre de bit
N = Nb / log2(M); % nombre de symbole
Db = 0.1;% débit binaire
D = Db / log2(M); % débit symbole
T = 1 / D; % période
NPTS = 1000; % nombre de points
NFFT = 2 ^nextpow2(NPTS);
f = linspace(-Fs/2,Fs/2,NFFT);
fc = 0.2; % fréquence de coupure

% Calcul
alpha_k=randint(Nb,1); % générations des bits aléatoirement
symbole=pskmod(alpha_k, M); % création des symboles aec pskmod
scatterplot(symbole); % affichage des symboles

%Signal Symboles
a = zeros(N * T, 1); % création d'un peignes de dirac avec un symbole toute les périodes T
for(i=0:N-1)
    a(i * T + 1) = symbole(i + 1);
end;

%Signal porte
porte=ones(T,1); 

%Signal Passe-Bas
a_filtre=filter(porte,1,a); % filtrage passe-bas à l'aide d'une porte

%Signal Passe-Bande
c = zeros(N * T, 1); 
for(i=1:N*T)
    c(i) = cos(2*pi*fc*i);
end;
s = a_filtre .* c; % transposition en fréquence à l'aide d'un cosinus de fréquence fc

% Affichage Temporel
subplot(3,1,1);
plot(real(s),'-ks','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',1);
hold on;
subplot(3,1,2);
plot(real(a_filtre),'-s','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',1);
subplot(3,1,3);
plot(imag(a_filtre),'-rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',1);


subplot(1, 1, 1);
% Affichage Densité Spectrale
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


%---------------------------------------------------------------------
% Partie Optionnelle



%---------------------------------------------------------------------