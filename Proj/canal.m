%------------------------
% Initialisation
%------------------------

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

Eg = T;
Va = 1;
Es = Eg*Va /2;
Eb = Es / log2(M);
No = Eb / (10^(3/20));

DSP_Bruit = zeros(1,nfft);
DSP_Demod_LP_Signal = zeros(1,nfft);
DSP_BPSignal = zeros(1,nfft);
DSP_Demod_BP_Signal_AWGN = zeros(1,nfft);

%-----------------------
%  PARTIE EMISSION
%-----------------------
for ii = 1 : 100

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

% Ajout du bruit blanc gaussien

Bruit = sqrt(No/2).*randn([1,(T/2)*Nb]);



%------------------------
% PARTIE RECEPTION
%------------------------

I_signal = BP_Signal.*real(Porteuse);
Q_signal = - BP_Signal.*imag(Porteuse);
I_signal_AWGN = (BP_Signal + Bruit).*real(Porteuse);
Q_signal_AWGN  = - (BP_Signal + Bruit).*imag(Porteuse);
Demod_BP_Signal = I_signal + j* Q_signal;       % Retour en bande de base
Demod_BP_Signal_AWGN  = I_signal_AWGN  + j* Q_signal_AWGN ; 


Demod_Porte = fliplr(Porte) / (T/2);               % Adaptation du filtre
Demod_LP_Signal = filter(Demod_Porte,1,Demod_BP_Signal);     % Filtrage
Demod_LP_Signal_AWGN  = filter(Demod_Porte,1,Demod_BP_Signal_AWGN);

Demod_Symboles = Demod_LP_Signal(T : T : N*T);
Demod_Symboles_AWGN = Demod_LP_Signal_AWGN(T : T : N*T);

Demod_Bits = pskdemod(Demod_Symboles, M);
Demod_Bits_AWGN  = pskdemod(Demod_Symboles_AWGN , M);
%Demod_Bits = pamdemod(Demod_Symboles, M);
%Demod_Bits = qamdemod(Demod_Symboles, M);
Symboles1 = de2bi(Demod_Bits,log2(M), Nb/log2(M)');



% Calcul des Densités Spectrales


DSP_Demod_BP_Signal_AWGN = DSP_Demod_BP_Signal_AWGN + fftshift((abs(fft(Demod_BP_Signal_AWGN,nfft)).^2) / nfft);            

DSP_LPSignal = fftshift((abs(fft(LP_Signal,nfft)).^2) / nfft);

DSP_BPSignal = DSP_BPSignal + fftshift((abs(fft(BP_Signal,nfft)).^2) / nfft);

DSP_Bruit = DSP_Bruit + fftshift((abs(fft(Bruit,nfft)).^2) / nfft);

DSP_Demod_LP_Signal = DSP_Demod_LP_Signal + fftshift((abs(fft(Demod_LP_Signal,nfft)).^2) / nfft); 

end

DSP_Bruit = DSP_Bruit / 100;
DSP_BPSignal = DSP_BPSignal /100;
DSP_Demod_LP_Signal = DSP_Demod_LP_Signal /100;
DSP_Demod_BP_Signal_AWGN = DSP_Demod_BP_Signal_AWGN /100;


%------------------------------
%  REPRESENTATION GRAPHIQUE
%------------------------------

scatterplot(Symboles)
scatterplot(Demod_Symboles)
scatterplot(Demod_Symboles_AWGN)

            % Représentation en temps

figure                          
Subplot(211)
plot(temps/T - 1/T, Demod_BP_Signal ,'r', temps/T - 1/T, Demod_BP_Signal_AWGN, 'b')
legend('Signal recu sans bruit','Signal recu avec bruit')
xlabel('[t/T]')
axis([0 N -20 20])

Subplot(212)
plot(temps/T - 1/T, Demod_LP_Signal ,'r', temps/T - 1/T, Demod_LP_Signal_AWGN, 'b')
legend('Signal recu sans bruit en sortie du filtre adapté','Signal recu avec bruit en sortie du filtre adapté')
xlabel('[t/T]')
axis([0 N -10 10])


            % Représentation en fréquence

Freq = linspace(-0.5,0.5,nfft);

figure                      
semilogy(Freq, DSP_BPSignal,'b', Freq, DSP_Demod_BP_Signal_AWGN, 'r', Freq, DSP_Bruit, 'g')
legend('Band-Pass Signal','Received Signal', 'AWGN')
xlabel('Normalized Frequency')
axis([-0.5 0.5 10^(-2) 10])