%-----------------------------------------
%  Initialisation                         
%-----------------------------------------
close all;
clear all;

%----------------------------------------
%Param�tres
%----------------------------------------

Fs = 1;
Nb = 10;     % Nb = 100 pour avoir la repr�sentation en fr�quence
M = 2;
T = 10;
N = Nb/log2(M);
D = 1/T;
Db = D * log2(M);
Tb = 1 / Db
Fc = 0.2;
temps = [1:N*T];
nfft = 2^nextpow2(N*T);

%-------------------------------------------
%  Modulation BPSK
%-------------------------------------------


Bits = randint(1,Nb);          % Generation d'une suite de Nb bits
Symboles = pskmod(Bits,M);     % Creation des symboles associ�s � cette suite

Peigne = zeros(1,N*T);         % Cr�ation du peigne
Peigne(1 : T : N*T) = Symboles(1 : 1 : N);
Porte = ones(1,T);             % Cr�ation de la porte

LP_Signal = filter(Porte,1,Peigne);   % Signal filtr� (= Signal Passe-Bas)

Porteuse = exp(j*2*pi*Fc*temps);      % Cr�ation de la porteuse

BP_Signal = real(LP_Signal.*Porteuse);% Transposition du signal autour de Fc

%-------------------------------------------
% Calcul des Densit�s Spectrales
%-------------------------------------------

DSP_Peigne = fftshift((abs(fft(Peigne,nfft)).^2) / nfft);

DSP_LPSignal = fftshift((abs(fft(LP_Signal,nfft)).^2) / nfft);

DSP_BPSignal = fftshift((abs(fft(BP_Signal,nfft)).^2) / nfft);

%--------------------------------------------
%  Affichage des courbes
%--------------------------------------------

scatterplot(Symboles)

% Repr�sentation en temps

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
plot(temps/T,imag(LP_Signal),'r', temps(1:1:N) - 1, imag(Peigne(1:T:N*T)),'go')
legend('Im(Low-Pass Signal)','Im(Symboles)')
xlabel('[t/T]')
axis([0 N -1 1])


% Repr�sentation en fr�quence

Freq = linspace(-0.5,0.5,nfft);

figure                      
semilogy(Freq, DSP_Peigne,'bo', Freq, DSP_LPSignal,'gd-', Freq, DSP_BPSignal, 'r-')
legend('Symboles','Low-Pass Signal','Band-Pass Signal')
xlabel('Normalized Frequency')
axis([-0.5 0.5 10^(-6) 100])
