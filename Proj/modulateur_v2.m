%  Initialisation

close all;
clear all;
Fs = 1;
Nb = 300;     % Nb = 100 pour avoir la représentation en fréquence
M = 8;
T = 10;
N = Nb/log2(M);
D = 1/T;
Db = D * log2(M);
Tb = 1 / Db
Fc = 0.2;
temps = [1:N*T];
nfft = 2^nextpow2(N*T);
phase_init = pi/4 % pi/4 pour le cas M=4

%  Modulation BPSK

Bits = randint(1,Nb);                       % Generation d'une suite de Nb bits
Symboles1 = bi2de(reshape(Bits,log2(M), Nb/log2(M))');
Symboles = pskmod(Symboles1,M,phase_init);                  % Creation des symboles associés à cette suite

Peigne = zeros(1,N*T);                                % Création d'un vecteur de 0
Peigne(1 : T : N*T) = Symboles(1 : 1 : N);            % Création du peigne

Porte = ones(1,T);                               % Création de la porte

LP_Signal = filter(Porte,1,Peigne);                % Signal filtré

Porteuse = exp(j*2*pi*Fc*temps);                    % Création de la porteuse

BP_Signal = real(LP_Signal.*Porteuse);             % Transposition autour de Fc

% Calcul des Densités Spectrales


DSP_Peigne = fftshift((abs(fft(Peigne,nfft)).^2) / nfft);               % Densités Spectrales

DSP_LPSignal = fftshift((abs(fft(LP_Signal,nfft)).^2) / nfft);

DSP_BPSignal = fftshift((abs(fft(BP_Signal,nfft)).^2) / nfft);


%  Affichage des courbes


scatterplot(Symboles)
axis([-1 1 -1 1])

            % Représentation en temps

figure                          
Subplot(311)
plot(temps/T, real(BP_Signal),'r')
legend('Band-Pass Signal')
xlabel('[t/T]')
Subplot(312)
plot(temps/T,real(LP_Signal), 'r',temps(1:1:N) - 1, real(Peigne(1:T:N*T)),'go')
legend('Re(Low-Pass Signal)','Re(Symboles)')
xlabel('[t/T]')
Subplot(313)
plot(temps/T,imag(LP_Signal),'r',temps(1:1:N) - 1, imag(Peigne(1:T:N*T)),'go')
legend('Im(Low-Pass Signal)','Im(Symboles)')
xlabel('[t/T]')
axis([0 N -1 1])


            % Représentation en fréquence

Freq = linspace(-0.5,0.5,nfft);

figure                      
semilogy(Freq, DSP_Peigne,'bo', Freq, DSP_LPSignal,'gd-', Freq, DSP_BPSignal, 'r-')
legend('Symboles','Low-Pass Signal','Band-Pass Signal')
xlabel('Normalized Frequency')
axis([-0.5 0.5 10^(-6) 100])