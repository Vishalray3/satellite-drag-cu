%% Change in a and x analytically
function [del_a,del_x, Cd_avg] = analytical_change_bff_iner(del_dash,a_sma,e,rho0,betaH,x0,An,order,case_ecc)


Cd_avg = 0;
bx = a_sma*e*betaH;
for nn = 1:order+1
    n = nn-1;
    if mod(n,2) == 0
        for k = 0:floor(n/2)
            Cd_avg = Cd_avg + An(nn)*(-1)^k*nchoosek(n,2*k)*(hm(n,k,0,bx) + ((k-n/2)*hm(n,k,0,bx)+n/2*hm(n,k,2,bx))*e^2);
        end
    else
        for k = 0:floor(n/2)
            Cd_avg = Cd_avg + An(nn)*(-1)^k*nchoosek(n,2*k)*(gm(n,k,0,bx) + ((k-n/2)*gm(n,k,0,bx)+n/2*gm(n,k,2,bx))*e^2);
        end
    end
end
Cd_avg = Cd_avg/(2*pi*besseli(0,bx));

switch case_ecc
    case 'low'
        bx = a_sma*e*betaH;
        int_term_a = 0;
        int_term_x = 0;
        for nn = 1:order+1
            n = nn-1;
            if mod(n,2) == 0
                for k = 0:floor(n/2)
                    int_term_a = int_term_a + An(nn)*(-1)^k*nchoosek(n,2*k)*(hm(n,k,0,bx) + 2*gm(n,k,1,bx)*e + ((k-n/2)*hm(n,k,0,bx)+1/2*(n+3)*hm(n,k,2,bx))*e^2 ...
                        + ((2*k-n)*gm(n,k,1,bx) + (n+1)*gm(n,k,3,bx))*e^3);
                    int_term_x = int_term_x + An(nn)*(-1)^k*nchoosek(n,2*k)*(gm(n,k,1,bx) + (hm(n,k,0,bx) + hm(n,k,2,bx))*e ...
                        + 1/2*((2+2*k-n)*gm(n,k,1,bx)+(n+1)*gm(n,k,3,bx))*e^2+ 1/2*((2*k-n)*hm(n,k,0,bx) + (2*k+1)*hm(n,k,2,bx) + (n+1)*hm(n,k,4,bx))*e^3);
                end
            else
                for k = 0:floor(n/2)
                    int_term_a = int_term_a + An(nn)*(-1)^k*nchoosek(n,2*k)*(gm(n,k,0,bx) + 2*hm(n,k,1,bx)*e + ((k-n/2)*gm(n,k,0,bx)+1/2*(n+3)*gm(n,k,2,bx))*e^2 ...
                        + ((2*k-n)*hm(n,k,1,bx) + (n+1)*hm(n,k,3,bx))*e^3);
                    int_term_x = int_term_x + An(nn)*(-1)^k*nchoosek(n,2*k)*(hm(n,k,1,bx) + (gm(n,k,0,bx) + gm(n,k,2,bx))*e ...
                        + 1/2*((2+2*k-n)*hm(n,k,1,bx)+(n+1)*hm(n,k,3,bx))*e^2+ 1/2*((2*k-n)*gm(n,k,0,bx) + (2*k+1)*gm(n,k,2,bx) + (n+1)*gm(n,k,4,bx))*e^3);
                end
            end
        end
        del_a =  -del_dash*a_sma^2*rho0*exp(betaH*(-x0))*int_term_a;
        del_x = -del_dash*a_sma^2*rho0*exp(betaH*(-x0))*int_term_x;
    case 'high'
        z = a_sma*e*betaH;
        int_term = 0;
        int_term2 = 0;
        for nn = 1:order+1
            n = nn - 1;
            nn_floor = floor(n/2)+1;
            for kk = 1:nn_floor
                k = kk - 1;
                L1 = ((-4*n + 6*k +1) - 8*e + e^2*(- 6*k + 3))/(4*(1 - e^2));
                L2 = (4*(2*n-3*k)^2+40*k-24*n+3 + 16*(4*n-6*k-1)*e - (72*k^2+32*k-24*n-48*k*n-50)*e^2 +16*(6*k+1)*e^3 ...
                    + (2*k-1)*(18*k+5)*e^4)/(32*(1-e^2)^2);
                int_term = int_term + (-1)^k*An(nn)*nchoosek(n,2*k)*(2/z/(1-e^2))^k*(gamma((2*k+1)/2) + L1/z*gamma((2*k+3)/2) + L2/z^2*gamma((2*k+5)/2));
                
                N1 = -((4*n - 6*k +3) + e^2*(6*k +1))/(4*(1 - e^2));
                N2 = (4*(2*n-3*k)^2 + 8*(n-k) - 5 + 32*e -2*(36*k^2+16*k-28*n-24*k*n+7)*e^2 + 32*e^3 ...
                    + (36*k^2+40*k+3)*e^4)/(32*(1-e^2)^2);
                int_term2 = int_term2 + (-1)^k*An(nn)*nchoosek(n,2*k)*(2/z/(1-e^2))^k*(gamma((2*k+1)/2) + N1/z*gamma((2*k+3)/2) + N2/z^2*gamma((2*k+5)/2));
            end
        end
        
        del_a = -del_dash*a_sma^2*rho0*sqrt(2/z)*(1 + e)^(3/2)/sqrt(1 - e)*int_term;
        del_x = -del_dash*a_sma^2*rho0*sqrt(2/z)*(1 + e)^(3/2)/sqrt(1 - e)*int_term2;
end


    function acoeff = a_coeff(k)
        acoeff = 1/2^(k)*nchoosek(k,k/2);
    end

    function h = hm(n,k,m,bx)
        p = n+m-2*k;
        S2 = 0;
        a_coeff1 = 0;
        if mod(p,2) == 0
            a_coeff1 = a_coeff(p);
            for i = 0:p/2-1
                S2 = S2 + 1/2^(p-1)*nchoosek(p,i)*besseli(p-2*i,bx);
            end
        end
        p = 2*k;
        S1 = 0;
        for j = 0:p/2-1
            S1 = S1 + (-1)^(p/2)/2^(p-1)*(-1)^j*nchoosek(p,j)*besseli(p-2*j,bx);
        end
        p1 = 2*k;
        p2 = n+m-2*k;
        S12 = 0;
        if mod(p2,2) == 0
            for j = 0:p1/2-1
                for i = 0:p2/2-1
                    S12 = S12 + (-1)^(p1/2)/2^(p1-1)/2^(p2-1)*(-1)^j*nchoosek(p1,j)*nchoosek(p2,i)*(besseli(p1+p2-2*i-2*j,bx) +besseli(p2-p1-2*i+2*j,bx));
                end
            end
        end
        
        h = 2*pi*(a_coeff(2*k)*a_coeff1*besseli(0,bx) + a_coeff(2*k)*S2 + a_coeff1*S1 + S12/2);
        
    end

    function g = gm(n,k,m,bx)
        p = n+m-2*k;
        S3 = 0;
        if mod(p,2) ~= 0
            for i = 0:(p-1)/2
                S3 = S3 + 1/2^(p-1)*nchoosek(p,i)*besseli(p-2*i,bx);
            end
        end
        p1 = 2*k;
        p2 = n+m-2*k;
        S13 = 0;
        if mod(p2,2) ~= 0
            for j = 0:p1/2-1
                for i = 0:(p2-1)/2
                    S13 = S13 + (-1)^(p1/2)/2^(p1-1)/2^(p2-1)*(-1)^j*nchoosek(p1,j)*nchoosek(p2,i)*(besseli(p1+p2-2*i-2*j,bx) +besseli(p2-p1-2*i+2*j,bx));
                end
            end
        end
        
        g = 2*pi*(a_coeff(2*k)*S3 + S13/2);
    end

end