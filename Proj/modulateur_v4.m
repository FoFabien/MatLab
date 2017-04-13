%  Initialisation

close all;
clear all;
Fs = 1;
Nb = 100;     % Nb = 100 pour la repr�sentation en fr�quence
M = 2;
T = 10;
N = Nb/log2(M);
D = 1/T;
Db = D * log2(M);
Tb = 1 / Db
Fc = 0.2;
temps = [1:N*T];
nfft = 2^nextpow2(N*T);
Nstat = 100;    % Nstat = 100 pour la repr�sentation en fr�quence
DSP_Signal_FSK_MATLAB = zeros(1,nfft);
DSP_Signal_FSK_MATLAB2 = zeros(1,nfft);
DSP_Signal_FSK = zeros(1,nfft);

%  Modulation BFSK avec la fonction fskmod (sans continuit� de phase)

for jj = 1 : Nstat
    
Bits = randint(1,Nb);                       % Generation d'une suite de Nb bits
symboles = bi2de(reshape(Bits,log2(M), Nb/log2(M))');
Symboles = fskmod(symboles,M,1/(2*T),T,1,'discont');                  % Creation des symboles associ�s � cette suite
%scatterplot(Symboles)
Peigne = zeros(1,N*T);                                % Cr�ation d'un vecteur de 0
Peigne(1 : T : N*T) = Symboles(1 : 1 : N);            % Cr�ation du peigne

Porte = ones(1,T);                               % Cr�ation de la porte

LP_Signal = filter(Porte,1,Peigne);                % Signal filtr�

Porteuse = exp(j*2*pi*Fc*temps);                    % Cr�ation de la porteuse

BP_Signal = real(LP_Signal.*Porteuse);             % Transposition autour de Fc

Signal_FSK_MATLAB = real(Porteuse.*Symboles');


% Modulation BFSK avec la fonction fskmod (avec continuit� de phase)

Bits2 = randint(1,Nb);                       % Generation d'une suite de Nb bits
symboles2 = bi2de(reshape(Bits2,log2(M), Nb/log2(M))');
Symboles2 = fskmod(symboles2,M,1/(2*T),T,1,'cont');                  % Creation des symboles associ�s � cette suite

Peigne2 = zeros(1,N*T);                                % Cr�ation d'un vecteur de 0
Peigne2(1 : T : N*T) = Symboles2(1 : 1 : N);            % Cr�ation du peigne

Porte2 = ones(1,T);                               % Cr�ation de la porte

LP_Signal2 = filter(Porte2,1,Peigne2);                % Signal filtr�

Porteuse2 = exp(j*2*pi*Fc*temps);                    % Cr�ation de la porteuse

BP_Signal2 = real(LP_Signal2.*Porteuse2);             % Transposition autour de Fc

Signal_FSK_MATLAB2 = real(Porteuse2.*Symboles2');



%  Modulation BFSK sans la fonction fskmod (sans continuit� de phase)

F1= cos(2*pi*(Fc-1/(4*T))*temps);
F2= cos(2*pi*(Fc+1/(4*T))*temps);
Signal_FSK = [];

for ii = 1 : N
    
    Signal_FSK =[Signal_FSK ((1-symboles(ii)).*F2((ii-1)*T + 1 : ii*T ) + symboles(ii).*F1((ii-1)*T +1 : ii*T))];
    
end


% Calcul des Densit�s Spectrales


DSP_Signal_FSK_MATLAB = DSP_Signal_FSK_MATLAB + fftshift((abs(fft(Signal_FSK_MATLAB,nfft)).^2) / nfft);

DSP_Signal_FSK_MATLAB2 = DSP_Signal_FSK_MATLAB2 + fftshift((abs(fft(Signal_FSK_MATLAB2,nfft)).^2) / nfft);

DSP_Signal_FSK = DSP_Signal_FSK + fftshift((abs(fft(Signal_FSK,nfft)).^2) / nfft);


end

DSP_Signal_FSK_MATLAB = DSP_Signal_FSK_MATLAB / Nstat;
DSP_Signal_FSK_MATLAB2 = DSP_Signal_FSK_MATLAB2 / Nstat;
DSP_Signal_FSK = DSP_Signal_FSK / Nstat;

%  Affichage des courbes pour la comparaison "Avec fskmod / Sans fskmod"

            % Repr�sentation en temps

figure                          
Subplot(211)
plot(temps/T, Signal_FSK ,'r')
legend('BFSK disc phase')
xlabel('[t/T]')

Subplot(212)
plot(temps/T, Signal_FSK_MATLAB, 'b')
legend('BFSK MATLAB disc phase')
xlabel('[t/T]')



            % Repr�sentation en fr�quence

Freq = linspace(-0.5,0.5,nfft);

figure                      
semilogy(Freq, DSP_Signal_FSK,'bo', Freq, DSP_Signal_FSK_MATLAB,'g-')
legend('BFSK disc phase','BFSK MATLAB disc phase')
xlabel('Normalized Frequency')
axis([-0.5 0.5 0.01 100])




%  Affichage des courbes pour la comparaison "Continuit� de phase / Discontinuit� de phase"

            % Repr�sentation en temps

figure                          
Subplot(211)
plot(temps/T, Signal_FSK_MATLAB2 ,'r')
legend('BFSK MATLAB cont phase')
xlabel('[t/T]')

Subplot(212)
plot(temps/T, Signal_FSK_MATLAB, 'b')
legend('BFSK MATLAB disc phase')
xlabel('[t/T]')



            % Repr�sentation en fr�quence

Freq = linspace(-0.5,0.5,nfft);

figure                      
semilogy(Freq, DSP_Signal_FSK_MATLAB2 ,'b-', Freq, DSP_Signal_FSK_MATLAB,'r-')
legend('BFSK MATLAB cont phase','BFSK MATLAB disc phase')
xlabel('Normalized Frequency')
axis([-0.5 0.5 0.0001 100])




