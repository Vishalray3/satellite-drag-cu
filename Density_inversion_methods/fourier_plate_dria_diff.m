function An = fourier_plate_dria_diff(Ai, Ar, s, r,Ci, order,flag_diff)

for jj = 1:numel(Ai)
    In = besseli([0:floor(order/2)+2],-s^2*Ci(jj)^2/2);
    In1 = besseli([-1:floor(order/2)+1],-s^2*Ci(jj)^2/2);
    In2 = besseli([1:floor(order/2)+3],-s^2*Ci(jj)^2/2);
    Ind = 1/2*(In1+In2)*(-s*Ci(jj)^2);
    for nn = 1:order+1
        n = nn-1;
        switch flag_diff
            case 's'
                if n == 0
                    An(nn,jj) = 0;
                elseif n == 1
                    An(nn,jj) = Ai(jj)/(pi*Ar)*(-1/s^3*Ci(jj)*pi + r(jj)*pi/12*Ci(jj)^3*((exp(-s^2*Ci(jj)^2/2)-s^2*Ci(jj)^2*exp(-s^2*Ci(jj)^2/2))...
                        *(9*In(1)-8*In(2)-In(3)) + s*exp(-s^2*Ci(jj)^2/2)*(9*Ind(1)-8*Ind(2)-Ind(3)))+ r(jj)*Ci(jj)*pi/(2*s)*exp(-s^2*Ci(jj)^2/2)*(Ind(2)+Ind(1))...
                        +r(jj)*Ci(jj)*pi/2*(In(2)+In(1))*(-exp(-s^2*Ci(jj)^2/2)/s^2 - Ci(jj)^2*exp(-s^2*Ci(jj)^2/2)));
                elseif n == 2
                    An(nn,jj) = Ai(jj)/(pi*Ar)*(2*sqrt(pi)/s*exp(-s^2*Ci(jj)^2/2)*Ind(2) + 2*sqrt(pi)/s*In(2)*(-1/s^2*exp(-s^2*Ci(jj)^2/2)...
                        -Ci(jj)^2*exp(-s^2*Ci(jj)^2/2)) + (-1/s^2*sqrt(pi)*Ci(jj)^2*exp(-s^2*Ci(jj)^2/2) + (1+1/(2*s^2))*sqrt(pi)*Ci(jj)^2*(exp(-s^2*Ci(jj)^2/2)...
                        -s^2*Ci(jj)^2*exp(-s^2*Ci(jj)^2/2)))*(3*In(1)-2*In(2)-In(3))/3 + ...
                        (1+1/(2*s^2))*sqrt(pi)*Ci(jj)^2*s*exp(-s^2*Ci(jj)^2/2)*(3*Ind(1)-2*Ind(2)-Ind(3))/3);
                elseif n > 2 && mod(n,2) == 0
                    k = n/2;
                    An(nn,jj) = Ai(jj)/(pi*Ar)*(2*sqrt(pi)/s*exp(-s^2*Ci(jj)^2/2)*Ind(k+1) + 2*sqrt(pi)*In(k+1)*(-1/s^2*exp(-s^2*Ci(jj)^2/2)...
                        -Ci(jj)^2*exp(-s^2*Ci(jj)^2/2))+ (1+1/(2*s^2))*sqrt(pi)*Ci(jj)^2*s*exp(-s^2*Ci(jj)^2/2)*((Ind(k+1) - Ind(k+2))/(2*k+1) ...
                        + (Ind(k) - Ind(k+1))/(2*k-1)) + sqrt(pi)*Ci(jj)^2*((In(k+1) - In(k+2))/(2*k+1) + (In(k) - In(k+1))/(2*k-1))...
                        *(exp(-s^2*Ci(jj)^2/2)/s^2 + (1+1/(2*s^2))*exp(-s^2*Ci(jj)^2/2)*(1-s^Ci(jj)^2)));
                elseif n > 2 && mod(n,2) ~= 0
                    k = (n-1)/2;
                    An(nn,jj) = Ai(jj)/(pi*Ar)*(r(jj)*pi/2*Ci(jj)^3*(exp(-s^2*Ci(jj)^2/2)-s^2*Ci(jj)^2*exp(-s^2*Ci(jj)^2/2))*((In(k+1) - In(k+2))/(2*k+1) ...
                        + (In(k+2) - In(k+3))/(4*k+6) + (In(k) - In(k+1))/(4*k-2)) + r(jj)*pi/2*Ci(jj)^3*s*exp(-s^2*Ci(jj)^2/2)*((Ind(k+1) - Ind(k+2))/(2*k+1) ...
                        + (Ind(k+2) - Ind(k+3))/(4*k+6) + (Ind(k) - Ind(k+1))/(4*k-2))+ r(jj)*Ci(jj)*pi/(2*s)*exp(-s^2*Ci(jj)^2/2)*(Ind(k+2) + Ind(k+1))...
                        + r(jj)*Ci(jj)*pi*(In(k+2) + In(k+1))*(-1/(2*s^2)*exp(-s^2*Ci(jj)^2/2)-Ci(jj)^2*exp(-s^2*Ci(jj)^2/2)));
                end
                
            case {'r_s', 'r_ads'}
                if n == 0
                    An(nn,jj) = 0;
                elseif n == 1
                    An(nn,jj) = Ai(jj)/(pi*Ar)*(pi/12*s*Ci(jj)^3*exp(-s^2*Ci(jj)^2/2)*(9*In(1)-8*In(2)-In(3)) + Ci(jj)*pi/(2*s)*exp(-s^2*Ci(jj)^2/2)*(In(2)+In(1)));
                elseif n == 2
                    An(nn,jj) = Ai(jj)/(pi*Ar)*(Ci(jj)^2*pi^(3/2)/4);
                elseif n > 2 && mod(n,2) == 0
                    k = n/2;
                    An(nn,jj) = 0;
                elseif n > 2 && mod(n,2) ~= 0
                    k = (n-1)/2;
                    An(nn,jj) = Ai(jj)/(pi*Ar)*(pi/2*s*Ci(jj)^3*exp(-s^2*Ci(jj)^2/2)*((In(k+1) - In(k+2))/(2*k+1)+ (In(k+2) - In(k+3))/(4*k+6)...
                        + (In(k) - In(k+1))/(4*k-2)) + Ci(jj)*pi/(2*s)*exp(-s^2*Ci(jj)^2/2)*(In(k+2) + In(k+1)));
                end 

        end
        
    end
end