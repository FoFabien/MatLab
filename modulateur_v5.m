function sd=modulateur_v5(M, continuite)

%Variable
Fs = 1;
Nb = 10; % 10 en temporel
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

% là encore, même chose que dans les modulateur_v4. On indique le M et la continuité en paramètre

% Calcul
alpha_k=randint(ceil(Nb / nextpow2(M)),nextpow2(M));
entiers = bi2de(alpha_k);

% modulation avec FSKMOD
sd=fskmod(entiers, M, delta, T, Fs,continuite);

e = zeros(N*T, 1);
for(ii=1:N*T)
   e(ii,1) =  exp(j * 2 * pi * fc * (ii - 1));
end;

plot(real(sd .* e));

semilogy(f, fftshift(abs(fft(real(sd .* e), 2^nextpow2(npts))).^2 /NFFT));