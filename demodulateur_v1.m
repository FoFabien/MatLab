function sd=demodulateur_v1(M)

%Variable
Fs = 1;
Nb = 100;
N = Nb / log2(M);
Db = 0.1;
D = Db / log2(M);
T = 1 / D;
npts = N*T;
NFFT = 2 ^nextpow2(npts);
f = linspace(-Fs/2,Fs/2,NFFT);
fc = 0.2;
delta = 1 / (2*T);
phi = 0;


signal_recu=modulateur(M, 'psk');

% création des deux porteuses pures et retour en bande de base
for(t1=0:N-1)
s1(t1+1,1)= signal_recu(t1+1,1) * cos(2 * pi * fc * t1 + phi);
s2(t1+1,1)= - signal_recu(t1+1,1) * sin(2 * pi * fc * t1 + phi);
end;


% filtrage adapté
adapt_porte = fliplr(ones(T,1)); % création du filtrage adapté au filtre de mise en forme
adapt1 = filter(adapt_porte, 1, s1); 
adapt2 = filter(adapt_porte, 1, s2);

% échantillonne au rythme symbole T
s_echantillon1 = adapt1(T:T:N,1);
s_echantillon2 = adapt2(T:T:N,1);

%s_final = s_echantillon1 + 1i * s_echantillon2;
s_final = adapt1 + 1i * adapt2;

s_demod = pskdemod(s_final, M);

sd = de2bi(s_demod, nextpow2(M)); % les bits obtenus en sortie ne correspondent pas à ceux en avant la modulation. On soupçonne une problème sur le filtrage

%représentation fréquentielle -------------------------------------------
figure
semilogy(f, fftshift(abs(fft(real(adapt1 + 1i * adapt2), 2^nextpow2(npts))).^2 /NFFT));
hold on
semilogy(f, fftshift(abs(fft(signal_recu, 2^nextpow2(npts))).^2 /NFFT), '-r');