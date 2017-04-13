%  Initialisation

close all;
clear all;
Fs = 1;
Nb = 10;     % Nb = 100 pour la représentation en fréquence
M = 4;
T = 10;
N = Nb/log2(M);
D = 1/T;
Db = D * log2(M);
Tb = 1 / Db
Fc = 0.2;
temps = [1:N*T];
nfft = 2^nextpow2(N*T);
Nstat = 1;    % Nstat = 100 pour la représentation en fréquence
DSP_Signal_FSK_MATLAB = zeros(1,nfft);
DSP_Signal_FSK_MATLAB2 = zeros(1,nfft);

%  Modulation BFSK avec la fonction fskmod (sans continuité de phase)

for jj = 1 : Nstat
    
Bits = randint(1,Nb);                       % Generation d'une suite de Nb bits
symboles = bi2de(reshape(Bits,log2(M), Nb/log2(M))');
Symboles = fskmod(symboles,M,1/(2*T),T,1,'discont');                  % Creation des symboles associés à cette suite

Peigne = zeros(1,N*T);                                % Création d'un vecteur de 0
Peigne(1 : T : N*T) = Symboles(1 : 1 : N);            % Création du peigne

Porte = ones(1,T);                               % Création de la porte

LP_Signal = filter(Porte,1,Peigne);                % Signal filtré

Porteuse = exp(j*2*pi*Fc*temps);                    % Création de la porteuse

BP_Signal = real(LP_Signal.*Porteuse);             % Transposition autour de Fc

Signal_FSK_MATLAB = real(Porteuse.*Symboles');


% Modulation BFSK avec la fonction fskmod (avec continuité de phase)

Bits2 = randint(1,Nb);                       % Generation d'une suite de Nb bits
symboles2 = bi2de(reshape(Bits2,log2(M), Nb/log2(M))');
Symboles2 = fskmod(symboles2,M,1/(2*T),T,1,'cont');                  % Creation des symboles associés à cette suite

Peigne2 = zeros(1,N*T);                                % Création d'un vecteur de 0
Peigne2(1 : T : N*T) = Symboles2(1 : 1 : N);            % Création du peigne

Porte2 = ones(1,T);                               % Création de la porte

LP_Signal2 = filter(Porte2,1,Peigne2);                % Signal filtré

Porteuse2 = exp(j*2*pi*Fc*temps);                    % Création de la porteuse

BP_Signal2 = real(LP_Signal2.*Porteuse2);             % Transposition autour de Fc

Signal_FSK_MATLAB2 = real(Porteuse2.*Symboles2');


% Calcul des Densités Spectrales


DSP_Signal_FSK_MATLAB = DSP_Signal_FSK_MATLAB + fftshift((abs(fft(Signal_FSK_MATLAB,nfft)).^2) / nfft);

DSP_Signal_FSK_MATLAB2 = DSP_Signal_FSK_MATLAB2 + fftshift((abs(fft(Signal_FSK_MATLAB2,nfft)).^2) / nfft);

end

DSP_Signal_FSK_MATLAB = DSP_Signal_FSK_MATLAB / Nstat;
DSP_Signal_FSK_MATLAB2 = DSP_Signal_FSK_MATLAB2 / Nstat;


%  Affichage des courbes pour la comparaison "Continuité de phase / Discontinuité de phase"

            % Représentation en temps

figure                          
Subplot(211)
plot(temps/T, Signal_FSK_MATLAB2 ,'r')
legend('BFSK MATLAB cont phase')
xlabel('[t/T]')

Subplot(212)
plot(temps/T, Signal_FSK_MATLAB, 'r')
legend('BFSK MATLAB disc phase')
xlabel('[t/T]')



            % Représentation en fréquence

Freq = linspace(-0.5,0.5,nfft);

figure                      
semilogy(Freq, DSP_Signal_FSK_MATLAB2 ,'b-', Freq, DSP_Signal_FSK_MATLAB,'r-')
legend('BFSK MATLAB cont phase','BFSK MATLAB disc phase')
xlabel('Normalized Frequency')
axis([-0.5 0.5 0.0001 100])




