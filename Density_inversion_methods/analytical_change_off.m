%% Change in a and x analytically
function [del_a,del_x] = analytical_change_off(del_dash,a_sma,e,rho0,betaH,x0,An,order,case_ecc)

switch case_ecc
    case 'low'
        bx = a_sma*e*betaH;
        for nn = 1:order+5
            In(nn,1) = besseli(nn-1, bx);
        end
        In0 = In(1:order+1);
        In1p = In(2:order+2);
        In2p = In(3:order+3);
        In3p = In(4:order+4);
        In4p = In(5:order+5);
        if order > 0
            In1n = [In(2);In(1:order)];
            
            In2n = [In(3);In(2);In(1:order-1)];
            if order > 1
                In3n = [In(4);In(3);In(2);In(1:order-2)];
            else
                In3n = [In(4);In(3)];
            end
            if order > 2
                In4n = [In(5);In(4);In(3);In(2);In(1:order-3)];
            elseif order > 1
                In4n = [In(5);In(4);In(3)];
            else
                In4n = [In(5);In(4)];
            end
        else
            In1n = In1p;
            In2n = In2p;
            In3n = In3p;
            In4n = In4p;
        end
        
        del_a = -2*pi*del_dash*a_sma^2*rho0*exp(betaH*(-x0))*(An*In0 + e*An*(In1p+In1n) + 3/8*e^2*An*(2*In0+In2p+In2n) + 1/8*e^3*An*(3*(In1p+In1n) + (In3p+In3n)));
        
        del_x = -2*pi*del_dash*a_sma^2*rho0*exp(betaH*(-x0))*(An/2*(In1p+In1n) + e/2*An*(3*In0 + 1/2*(In2p+In2n)) + e^2/16*An*(11*(In1p+In1n)+(In3p+In3n))...
            + e^3/16*An*(7*In0 + 4*(In2p+In2n) + 1/2*(In4p+In4n)));
    case 'high'
        z = a_sma*e*betaH;
        int_term = 0;
        int_term2 = 0;
        for nn = 1:order+1
            n = nn - 1;
            nn_floor = floor(n/2)+1;
            for kk = 1:nn_floor
                k = kk - 1;
                K1 = ((-4*n + 6*k +1) - 8*e + e^2*(4*n - 6*k + 3))/(4*(1 - e^2));
                K2 = ((4*n-6*k)*(4*n-6*k-6) + (4*k+3)+ 16*(4*n-6*k-1)*e + ((4*n-6*k)*(-8*n+12*k+4)-(8*k-50))*e^2 -16*(4*n-6*k-1)*e^3 ...
                       + ((4*n-6*k)*(4*n-6*k+2) +(4*k-5))*e^4)/(32*(1-e^2)^2);
                int_term = int_term + (-1)^k*An(nn)*nchoosek(n,2*k)*(2/z)^k*(gamma((2*k+1)/2) + K1/z*gamma((2*k+3)/2) + K2/z^2*gamma((2*k+5)/2));
                
                M1 = ((-4*n + 6*k -3) + e^2*(4*n - 6*k -1))/(4*(1 - e^2));
                M2 = ((4*n-6*k)*(4*n-6*k+2) + (4*k-5)+ 32*e -2*((4*n-6*k)*(4*n-6*k-2)+(4*k+7))*e^2 + 32*e^3 ...
                       + ((4*n-6*k)*(4*n-6*k-6) +(4*k+3))*e^4)/(32*(1-e^2)^2);  
                int_term2 = int_term2 + (-1)^k*An(nn)*nchoosek(n,2*k)*(2/z)^k*(gamma((2*k+1)/2) + M1/z*gamma((2*k+3)/2) + M2/z^2*gamma((2*k+5)/2));
            end
        end
        
        del_a = -del_dash*a_sma^2*rho0*sqrt(2/z)*(1 + e)^(3/2)/sqrt(1 - e)*int_term;
        del_x = -del_dash*a_sma^2*rho0*sqrt(2/z)*(1 + e)^(3/2)/sqrt(1 - e)*int_term2;
end
end