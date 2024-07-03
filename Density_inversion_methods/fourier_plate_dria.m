function An = fourier_plate_dria(Ai, Ar, s, r,Ci, order)

for jj = 1:numel(Ai)
In = besseli([0:floor(order/2)+2],-s^2*Ci(jj)^2/2);
for nn = 1:order+1
    n = nn-1;
    if n == 0
        An(nn,jj) = Ai(jj)/(2*pi*Ar)*(2*sqrt(pi)/s*exp(-s^2*Ci(jj)^2/2)*In(1) + 2*(1+1/(2*s^2))*sqrt(pi)*Ci(jj)^2*s*exp(-s^2*Ci(jj)^2/2)*(In(1) - In(2)) ...
            +r(jj)*Ci(jj)^2*pi^(3/2)/2);
    elseif n == 1
        An(nn,jj) = Ai(jj)/(pi*Ar)*((1+1/(2*s^2))*Ci(jj)*pi + r(jj)*pi/12*s*Ci(jj)^3*exp(-s^2*Ci(jj)^2/2)*(9*In(1)-8*In(2)-In(3)) ...
            + r(jj)*Ci(jj)*pi/(2*s)*exp(-s^2*Ci(jj)^2/2)*(In(2)+In(1)));
    elseif n == 2
        An(nn,jj) = Ai(jj)/(pi*Ar)*(2*sqrt(pi)/s*exp(-s^2*Ci(jj)^2/2)*In(2) + (1+1/(2*s^2))*sqrt(pi)*Ci(jj)^2*s*exp(-s^2*Ci(jj)^2/2)*(3*In(1)-2*In(2)-In(3))/3 ...
            +r(jj)*Ci(jj)^2*pi^(3/2)/4);
    elseif n > 2 && mod(n,2) == 0
        k = n/2;
        An(nn,jj) = Ai(jj)/(pi*Ar)*(2*sqrt(pi)/s*exp(-s^2*Ci(jj)^2/2)*In(k+1) + (1+1/(2*s^2))*sqrt(pi)*Ci(jj)^2*s*exp(-s^2*Ci(jj)^2/2)*((In(k+1) - In(k+2))/(2*k+1) ...
             + (In(k) - In(k+1))/(2*k-1)));
    elseif n > 2 && mod(n,2) ~= 0
        k = (n-1)/2;
        An(nn,jj) = Ai(jj)/(pi*Ar)*(r(jj)*pi/2*s*Ci(jj)^3*exp(-s^2*Ci(jj)^2/2)*((In(k+1) - In(k+2))/(2*k+1) + (In(k+2) - In(k+3))/(4*k+6) + (In(k) - In(k+1))/(4*k-2))...
            + r(jj)*Ci(jj)*pi/(2*s)*exp(-s^2*Ci(jj)^2/2)*(In(k+2) + In(k+1)));
    end  
end
end