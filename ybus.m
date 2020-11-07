% Programme pour l'admittance du bus

function Y = ybus(num)  % Return Y

linedata = linedatas( );        % la fonction du dataline
fb = linedata(:,1);             % From bus number...
tb = linedata(:,2);             % To bus number...
r = linedata(:,3);              % Resistance, R...
x = linedata(:,4);              % Reactance, X...
b = linedata(:,5);              % Ground Admittance, B/2...
a = linedata(:,6);              % Tap setting value..
z = r + i*x;                    % z matrice...
y = 1./z;                       % admittance de chaque element
b = i*b;                        

nb = max(max(fb),max(tb));      % No. of buses...
nl = length(fb);                % No. of branches...
Y = zeros(nb,nb);               
 
 % les elements qui n'appartiennt pas au diagonal
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

 