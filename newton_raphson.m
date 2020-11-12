function [Volt,theta] = newton_raphson(Y,busd)
nbus=14;
bus = busd(:,1);            % Bus Number..
type = busd(:,2);           % le Type de Bus 1-Slack, 2-PV, 3-PQ..
V = busd(:,3);              % la tension de chaque bus (pour les bus de type 3 v=1 pu comme valeur initiale)
del = busd(:,4);            % la phase de chaque tension
Pg = busd(:,5);        % PGi..
Qg = busd(:,6);        % QGi..
Pl = busd(:,7);        % PLi..
Ql = busd(:,8);        % QLi..
Qmin = busd(:,9);      % Minimum Reactive Power Limit..
Qmax = busd(:,10);     % Maximum Reactive Power Limit..
P = Pg - Pl;                % Pi = PGi - PLi..
Q = Qg - Ql;                % Qi = QGi - QLi..
Psp = P;                    % P Specified..
Qsp = Q;                    % Q Specified..
G = real(Y);                % Y=G+iB
B = imag(Y);              

pv = find(type == 2 | type == 1);   % PV Buses..
pq = find(type == 3);               % PQ Buses..
npv = length(pv);                   % No. of PV buses..
npq = length(pq);                   % No. of PQ buses..

Tol = 1;  
Iter = 1;
while (Tol > 1e-5)   
    
    P = zeros(14,1);
    Q = zeros(14,1);
    % determiner P et Q
    for i = 1:14
        for k = 1:14
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end

    % Vérification des violations des limites de Q
    if Iter <= 7 && Iter > 2    
        for n = 2:14
            if type(n) == 2
                QG = Q(n)+Ql(n);
                if QG < Qmin(n)
                    V(n) = V(n) + 0.01;
                elseif QG > Qmax(n)
                    V(n) = V(n) - 0.01;
                end
            end
         end
    end
    
    % DP et DQ
    dPa = Psp-P;
    dQa = Qsp-Q;
    k = 1;
    dQ = zeros(npq,1);
    for i = 1:14
        if type(i) == 3
            dQ(k,1) = dQa(i);
            k = k+1;
        end
    end
    dP = dPa(2:14);
    M = [dP; dQ];       % M=J*X
    
    % Jacobian
    % J  =    "  J11  J12  "
    %         "  J21  J22  "
    %
    %
    %
    % J1 - dPi/dtheta
    J1 = zeros(13,13);
    for i = 1:(13)
        m = i+1;
        for k = 1:(13)
            n = k+1;
            if n == m
                for n = 1:14
                    J1(i,k) = J1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
            else
                J1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    % J2 - Pi/dVi
    J2 = zeros(13,npq);
    for i = 1:(13)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k) = J2(i,k) + V(m)*G(m,m);
            else
                J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J3 - dQi/dtheta
    J3 = zeros(npq,13);
    for i = 1:npq
        m = pq(i);
        for k = 1:(13)
            n = k+1;
            if n == m
                for n = 1:14
                    J3(i,k) = J3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
            else
                J3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J4 - dQi/dVi
    J4 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:14
                    J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,k) = J4(i,k) - V(m)*B(m,m);
            else
                J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    J = [J1 J2; J3 J4]     % Jacobian Matrix..
    
    X = inv(J)*M;           % X=(J)^-1*M
    dTh = X(1:13);      % X=[dth dV]
    dV = X(14:end);       
    
    % Updating State Vectors..
    del(2:14) = dTh + del(2:14);    % Voltage Angle..
    k = 1;
    for i = 2:14
        if type(i) == 3
            V(i) = dV(k) + V(i);        % l'amplitude des bus de type 3
            k = k+1;
        end
    end
    
    Iter = Iter + 1;
    Tol = max(abs(M));                  % Tolerance..
    Volt=V;
    theta=del;
end


 subplot(211)
      bar(V);
      title('bus voltage magnitude')
      xlabel('bus number')
      ylabel('voltage in pu')
      
 subplot(212)
      bar((180/pi)*del);
     title('bus voltage phase')
     xlabel('bus number')
     ylabel('phase in degrees')
