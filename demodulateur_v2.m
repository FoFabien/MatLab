function sd=demodulateur_v2(M, type)

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

%----------------------------------------------------------------------------------------------------
%
%	Non test� � cause des soucis rencontr�s sur son pr�decesseur, demodulateur_v1.m
%
%----------------------------------------------------------------------------------------------------

signal_recu=modulateur(M, type);

% cr�ation des deux porteuses pures et retour en bande de base
for(t1=0:N-1)
s1(t1+1,1)= signal_recu(t1+1,1) * cos(2 * pi * fc * t1 + phi);
s2(t1+1,1)= - signal_recu(t1+1,1) * sin(2 * pi * fc * t1 + phi);
end;


% filtrage adapt�
adapt_porte = fliplr(ones(T,1)); % cr�ation du filtrage adapt� au filtre de mise en forme
adapt1 = filter(adapt_porte, 1, s1); 
adapt2 = filter(adapt_porte, 1, s2);

% �chantillonne au rythme symbole T
s_echantillon1 = adapt1(T:T:N,1);
s_echantillon2 = adapt2(T:T:N,1);

%s_final = s_echantillon1 + 1i * s_echantillon2;
s_final = adapt1 + 1i * adapt2;

if(type == 'psk') % choix de la d�modulation
s_demod = pskdemod(s_final, M);
elseif(type == 'pam')
s_demod = pamdemod(s_final, M);
elseif(type == 'qam')
s_demod = qamdemod(s_final, M);
elseif(type == 'fsk')
s_demod = fskdemod(s_final, M, delta, T, Fs);
end;

sd = de2bi(s_demod, nextpow2(M));

%repr�sentation fr�quentielle -------------------------------------------
figure
semilogy(f, fftshift(abs(fft(real(adapt1 + 1i * adapt2), 2^nextpow2(npts))).^2 /NFFT));
hold on
semilogy(f, fftshift(abs(fft(signal_recu, 2^nextpow2(npts))).^2 /NFFT), '-r');