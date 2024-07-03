function Cd_total = sphere_dria_cd(S, m_r, Ri, Tw, Vi, Alpha,frac,Kl,P_o,shape_model)
Vratio = sqrt(0.5*(1+Alpha*(4*Ri*Tw/Vi^2 -1)));
Cd_ads = (4*S^4 + 4*S^2 - 1)/(2*S^4)*erf(S) + (2*S^2+1)/(sqrt(pi)*S^3)*exp(-S^2) + 2*sqrt(pi)/3*Vratio;
Alpha_scham = 2.4*m_r/(1+m_r)^2;
Vratio = sqrt(0.5*(1+Alpha_scham*(4*Ri*Tw/Vi^2 -1)));
Cd_s = (4*S^4 + 4*S^2 - 1)/(2*S^4)*erf(S) + (2*S^2+1)/(sqrt(pi)*S^3)*exp(-S^2) + 2*sqrt(pi)/3*Vratio;
switch shape_model
    case 'sphere'
        Cd_total = frac*Cd_ads + (1-frac)*Cd_s;
    case 'sphere_jac'
        Cd_total = P_o/(1+Kl*P_o)^2*(Cd_ads - Cd_s);
end
