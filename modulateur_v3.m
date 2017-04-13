function sd=modulateur_v3(M)

%Variable
Fs = 1;
M = 16;
Nb = 1000;
N = Nb / log2(M);
Db = 0.1;
D = Db / log2(M);
T = 1 / D;
NPTS = 1000;
NFFT = 2 ^nextpow2(NPTS);
f = linspace(-Fs/2,Fs/2,NFFT);
fc = 0.2;

% Calcul
alpha_k=randint(ceil(Nb / nextpow2(M)),nextpow2(M)); % même chose que dans modulateur_v2.m
entiers = bi2de(alpha_k);

symbole_psk=pskmod(entiers, M); % cette fois on affiche toute les modulations différentes
figure
scatterplot(symbole_psk);

symbole_pam=pammod(entiers, M);
figure
scatterplot(symbole_pam);

symbole_qam=qammod(entiers, M);
figure
scatterplot(symbole_qam);