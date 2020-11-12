14%% Newton Raphson Load Flow analysis
clear all; close all; clc;


    nbus = 14;
    

Y = admittance_matrix( )             % appeler la matrice d'admittance
busd = busdatas( );      % les donnees de bus 14



disp('les resultats obtener par la methode NR ');

[voltage , theta] =newton_raphson(Y,busd)
disp(voltage);
disp(theta);
