%  Initialisation

close all;
clear all;
Fs = 1;
Nb = 10;     % Nb = 1000 pour avoir la représentation en fréquence
M = 2;
T = 10
N = Nb/log2(M)
D = 1/T;
Db = D * log2(M);
Tb = 1 / Db;
Fc = 0.2;
temps = [1:N*T];
nfft = 2^nextpow2(N*T);
phase_init = pi/4;

%  PARTIE EMISSION

Bits = randint(1,Nb);                       % Generation d'une suite de Nb bits
Symboles1 = bi2de(reshape(Bits,log2(M), Nb/log2(M))');
Symboles = fskmod(Symboles1,M,1/(2*T),T,1,'discont');                   % Creation des symboles associés à cette suit
Porteuse = exp(j*2*pi*Fc*(temps));

Signal_FSK_MATLAB = real(Porteuse.*Symboles');


% PARTIE RECEPTION

        % Avec fskdemod
        
env_compl = hilbert(Signal_FSK_MATLAB).*conj(Porteuse);
Demod_Bits = fskdemod(env_compl, M, 1/(2*T),T);

        % Sans fskdemod
       
F1= cos(2*pi*(Fc-1/(4*T))*temps);
F2= cos(2*pi*(Fc+1/(4*T))*temps);   
sortie0;
sortie1;
Signal_FSK1 = [];
Signal_FSK2 = [];

for ii = 1 : N
    
    rk1 = []
    Signal_FSK1 = [Signal_FSK1 (Signal_FSK_MATLAB.*F1((ii-1)*T + 1 : ii*T ))];
    Signal_FSK2 = [Signal_FSK2 (Signal_FSK_MATLAB.*F2((ii-1)*T + 1 : ii*T ))];
    
    Demod_Bits2 = max(Signal_FSK1,Signal_FSK2);
    
end
        
        
        
%  Affichage des courbes

scatterplot(Demod_Bits)
scatterplot(Demod_Bits2)

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
