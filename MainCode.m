14%% Newton Raphson Load Flow analysis
clear all; close all; clc;


    nbus = 14;
    

Y = admittance_matrix( )             % appeler la matrice d'admittance
busd = busdatas( );      % les donnees de bus 14



disp('les resultats obtener par la methode NR ');

[voltage , theta] =newton_raphson(Y,busd)
[Pi Qi Pg Qg Pl Ql]=loadflow(voltage,theta);
disp('Voltage:');
disp(voltage);
disp('angle: ');
disp(theta);
disp('elem active power');

disp(Pi);
disp('elem reactive power ');

disp(Qi);
disp('bus generator active power ');
disp(Pg);
disp('bus generator reactive power ');
disp(Qg);
disp('load active power ');
disp(Pl);
disp('load reactive power');
disp(Ql);
