%  Initialisation

close all;
clear all;
Fs = 1;
Nb = 10;     % Nb = 100 pour la repr�sentation en fr�quence
M = 4;
T = 10;
N = Nb/log2(M);
D = 1/T;
Db = D * log2(M);
Tb = 1 / Db
Fc = 0.2;
temps = [1:N*T];
nfft = 2^nextpow2(N*T);
Nstat = 1;    % Nstat = 100 pour la repr�sentation en fr�quence
DSP_Signal = zeros(1,nfft);

%  Modulation BFSK avec la fonction fskmod (sans continuit� de phase)

for jj = 1 : Nstat
    
% Modulation BFSK avec la fonction fskmod (avec continuit� de phase)

Bits = randint(1,Nb);                       % Generation d'une suite de Nb bits
symboles = bi2de(reshape(Bits,log2(M), Nb/log2(M))');
%Symboles = pammod(symboles,M);              % Choix de la modulation
%Symboles = qammod(symboles,M);
%Symboles = pskmod(symboles,M);
Symboles = fskmod(symboles,M,1/(2*T),T,1,'cont');                  % Creation des symboles associ�s � cette suite

Peigne = zeros(1,N*T);                                % Cr�ation d'un vecteur de 0
Peigne(1 : T : N*T) = Symboles(1 : 1 : N);            % Cr�ation du peigne

Porte = ones(1,T);                               % Cr�ation de la porte

LP_Signal = filter(Porte,1,Peigne);                % Signal filtr�

Porteuse = exp(j*2*pi*Fc*temps);                    % Cr�ation de la porteuse

BP_Signal = real(LP_Signal.*Porteuse);             % Transposition autour de Fc

Signal= real(Porteuse.*Symboles');


% Calcul des Densit�s Spectrales


DSP_Signal = DSP_Signal + fftshift((abs(fft(Signal,nfft)).^2) / nfft);

end

DSP_Signal = DSP_Signal / Nstat;


%  Affichage des courbes pour la comparaison "Continuit� de phase / Discontinuit� de phase"

            % Repr�sentation en temps

figure
plot(temps/T, Signal, 'r')
legend('Modulation choisie')
xlabel('[t/T]')



            % Repr�sentation en fr�quence

Freq = linspace(-0.5,0.5,nfft);

figure                      
semilogy(Freq, DSP_Signal ,'b-')
legend('Modulation choisie')
xlabel('Normalized Frequency')
axis([-0.5 0.5 0.0001 100])




