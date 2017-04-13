%  Initialisation

close all;
clear all;
Fs = 1;
Nb = 100;     % Nb = 100 pour avoir la représentation en fréquence
M = 2;
T = 10;
N = Nb/log2(M);
D = 1/T;
Db = D * log2(M);
Tb = 1 / Db;
Fc = 0.2;
temps = [1:N*T];
nfft = 2^nextpow2(N*T);
Nstat = 1;  % Nstat = 100 pour une meilleure estimation
DSP_Peigne = zeros(1,nfft);          % Définition des vecteurs Densités Spectrales
DSP_LPSignal = zeros(1,nfft);
DSP_BPSignal = zeros(1,nfft);

%  Modulation BPSK

for ii = 1 : Nstat
    
Bits = randint(1,Nb);                       % Generation d'une suite de Nb bits
Symboles = pskmod(Bits,M);                  % Creation des symboles associés à cette suite

Peigne = zeros(1,N*T);                                % Création d'un vecteur de 0
Peigne(1 : T : N*T) = Symboles(1 : 1 : N);            % Création du peigne

Porte = ones(1,T);                               % Création de la porte

LP_Signal = filter(Porte,1,Peigne);                % Signal filtré

Porteuse = exp(j*2*pi*Fc*temps);                    % Création de la porteuse

BP_Signal = real(LP_Signal.*Porteuse);             % Transposition autour de Fc

% Calcul des Densités Spectrales


DSP_Peigne = DSP_Peigne + fftshift((abs(fft(Peigne,nfft)).^2) / nfft);               % Densités Spectrales

DSP_LPSignal = DSP_LPSignal + fftshift((abs(fft(LP_Signal,nfft)).^2) / nfft);

DSP_BPSignal = DSP_BPSignal + fftshift((abs(fft(BP_Signal,nfft)).^2) / nfft);

end

DSP_Peigne = DSP_Peigne/Nstat;
DSP_LPSignal = DSP_LPSignal/Nstat;
DSP_BPSignal = DSP_BPSignal/Nstat;

% Courbes théoriques
Freq = linspace(-0.5,0.5,nfft);

DSP_Peigne_Theo = 1/T;
DSP_LPSignal_Theo = fftshift(abs(fft(T*Porte,nfft).^2)/nfft);
DSP_BPSignal_Theo = (T/4) * ((sinc((Freq - Fc)*T).^2) + (sinc((Freq + Fc)*T).^2));




%  Affichage des courbes

            % Représentation en fréquence




figure % Représentation en fréquence
Subplot(131)
semilogy(Freq, DSP_Peigne,'bo', Freq, DSP_Peigne_Theo,'r')
xlabel('Normalized Frequency')
legend('BPSK Symbols','PSD Symbols Theo')
axis([-0.5 0.5 10^(-6) 100])
Subplot(132)
semilogy(Freq, DSP_LPSignal,'bo', Freq, DSP_LPSignal_Theo,'r')
legend('Low-Pass Signal','PSD Low-Pass Theo')
xlabel('Normalized Frequency')
axis([-0.5 0.5 10^(-6) 100])
Subplot(133)
semilogy(Freq, DSP_BPSignal, 'bo', Freq, DSP_BPSignal_Theo,'r')
legend('Band-Pass Signal','PSD Band-Pass Theo')
xlabel('Normalized Frequency')
axis([-0.5 0.5 10^(-6) 100])