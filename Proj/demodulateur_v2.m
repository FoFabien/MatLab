%  Initialisation

close all;
clear all;
Fs = 1;
Nb = 500;     % Nb = 1000 pour avoir la représentation en fréquence
M = 4;
T = 10
N = Nb/log2(M)
D = 1/T;
Db = D * log2(M);
Tb = 1 / Db;
Fc = 0.2;
temps = [1:N*T];
nfft = 2^nextpow2(N*T);
phase_init = pi/4;

%-----------------------
%  PARTIE EMISSION
%-----------------------


Bits = randint(1,Nb);                       % Generation d'une suite de Nb bits
Symboles1 = bi2de(reshape(Bits,log2(M), Nb/log2(M))');
Symboles = pskmod(Symboles1,M,phase_init);                  % Creation des symboles associés à cette suite
%Symboles = pammod(Symboles1,M);
%Symboles = qammod(Symboles1,M);

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


%------------------------
% PARTIE RECEPTION
%------------------------


I_signal = BP_Signal.*real(Porteuse);
Q_signal = - BP_Signal.*imag(Porteuse);
Demod_BP_Signal = I_signal + j* Q_signal;       % Retour en bande de base

Demod_Porte = fliplr(Porte) / (T/2);               % Adaptation du filtre
Demod_LP_Signal = filter(Demod_Porte,1,Demod_BP_Signal);     % Filtrage
Demod_Symboles = Demod_LP_Signal(T : T : N*T);
Demod_Bits = pskdemod(Demod_Symboles, M);
%Demod_Bits = pamdemod(Demod_Symboles, M);
%Demod_Bits = qamdemod(Demod_Symboles, M);
Symboles1 = de2bi(Demod_Bits,log2(M), Nb/log2(M)');


DSP_Demod_LP_Signal = fftshift((abs(fft(Demod_LP_Signal,nfft)).^2) / nfft); 


%  Affichage des courbes

scatterplot(Symboles)
scatterplot(Demod_Symboles)

            % Représentation en temps

figure                          
Subplot(211)
plot(temps(1:1:N) - 1, real(Peigne(1:T:N*T)) ,'ro', temps/T - 1/T, real(LP_Signal), 'r-', temps(1:1:N) - 1/T, real(Demod_Symboles(1:1:N)), 'bo', temps/T - 1/T, real(Demod_LP_Signal),'b-')
legend('Re(Symbols)','Re(Low Pass Signal)', 'Re(Output Samples)', 'Re(Matched Filter Output)')
xlabel('[t/T]')
axis([0 N -20 20])
Subplot(212)
plot(temps(1:1:N) - 1, imag(Peigne(1:T:N*T)) ,'ro', temps/T - 1/T, imag(LP_Signal), 'r-', temps(1:1:N) - 1/T, imag(Demod_Symboles(1:1:N)), 'bo', temps/T - 1/T, imag(Demod_LP_Signal),'b-')
legend('Im(Symbols)','Im(Low Pass Signal)', 'Im(Output Samples)', 'Im(Matched Filter Output)')
xlabel('[t/T]')
axis([0 N -10 10])


            % Représentation en fréquence

Freq = linspace(-0.5,0.5,nfft);

figure                      
semilogy(Freq, DSP_BPSignal,'b', Freq, DSP_Demod_LP_Signal, 'r')
legend('Band Pass Signal','Matched Filter Output')
xlabel('Normalized Frequency')
axis([-0.5 0.5 10^(-6) 1000])