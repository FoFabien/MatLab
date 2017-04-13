function sd=modulateur(M, type) % là aussi, presque identique à modulateur_v5, on indique le type de modulation voulue en paramètre

%Variable
Fs = 1;
Nb = 100;
N = Nb / log2(M);
Db = 0.1;
D = Db / log2(M);
T = 1 / D;
NPTS = N*T;
NFFT = 2 ^nextpow2(NPTS);
f = linspace(-Fs/2,Fs/2,NFFT);
fc = 0.2;
delta = 1 / (2*T);

% Calcul
alpha_k = randint(ceil(Nb / nextpow2(M)),nextpow2(M))
entiers = bi2de(alpha_k);

if(type == 'psk')
s_mod=pskmod(entiers, M);
scatterplot(s_mod);
elseif(type == 'pam')
s_mod=pammod(entiers, M);
scatterplot(s_mod);
elseif(type == 'qam')
s_mod=qammod(entiers, M);
scatterplot(s_mod);
elseif(type == 'fsk')
s_mod=fskmod(entiers, M, delta, T, Fs);
scatterplot(s_mod);
end;

e = zeros(N, 1);
for(ii=1:N)
   e(ii,1) =  exp(1i * 2 * pi * fc * (ii - 1));
end;

figure
semilogy(f, fftshift(abs(fft(real(s_mod .* e), 2^nextpow2(NPTS))).^2 /NFFT)); % affichage de la représentation fréquentielle

%filtrage
porte=ones(T,1);
s_lowpass=filter(porte,T,s_mod); % signal passe-bas
for(t1=0:N-1)
sd(t1+1,1) = s_lowpass(t1+1,1) .* cos(2 * pi * fc * t1); % signal passe-bande
end;
