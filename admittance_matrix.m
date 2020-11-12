


            %_______________________________________________
            %                                               %
            % programme pour former la matrice d'admittances%
            %_______________________________________________%

function Y = ybus( )  % Returns Y

linedata = linedatas( );      % le tableau des donnees
fb = linedata(:,1);            % for bus
tb = linedata(:,2);             % to bus
r = linedata(:,3);              % Resistance, R
x = linedata(:,4);              % Reactance, X
b = linedata(:,5);              % l'admittance c
a = linedata(:,6);              
z = r + i*x;                    %la  matrice z
y = 1./z;                       % To get inverse of each element...
b = i*b;                        % Make B imaginary...

nb = 14      % No. of buses...
nl = length(fb);                % No. des  branches.
Y = zeros(nb,nb);               % Initialiser Y par des zeros
 
 % les elements qui n'apartiennent pas au diagonal
 for k = 1:nl
     Y(fb(k),tb(k)) =  - y(k)/a(k);
     Y(tb(k),fb(k)) = Y(fb(k),tb(k));
 end
 
 % les elements du diagonal
 for m = 1:nb
     for n = 1:nl
         if fb(n) == m
             Y(m,m) = Y(m,m) + y(n)/(a(n)^2) + b(n);
         elseif tb(n) == m
             Y(m,m) = Y(m,m) + y(n) + b(n);
         end
     end
 end
 