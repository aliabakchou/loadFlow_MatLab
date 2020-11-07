%% Newton Raphson Load Flow analysis
clear all; close all; clc;


    disp('using the IEEE 14 bus power system')

    nbus = 14;
    

Y = ybus(nbus) ;            
busd = busdatas( );      
BMva = 1;               



disp(Y);
