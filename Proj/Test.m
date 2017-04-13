%------------------------
% Initialisation
%------------------------

close all;
clear all;

%-------------
% Test
%-------------

N = 1000;
Y = randn([1,N]);
Z = zeros(1,N);
Auto_corr = conv(conj(Y),fliplr(Y))/N;

%----------------------------
% Représentation Graphique
%----------------------------

figure
plot([-N/2:1:N/2 - 1], Y, 'r')

figure
plot([-N+1:1:N-1], Auto_corr, 'b')
xlabel('Tau')
title('Auto-Correlation Functiun of an AWGN')


%-------------------
% Histogramme
%-------------------

M = 3;

for k = 1 : M
    
    Z = Z + rand([1,N]);
    
end

Auto_corr2 = conv(conj(Z),fliplr(Z))/N;

figure
hist(Z,50);
title('Histogram of an AWGN')