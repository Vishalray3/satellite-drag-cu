function [Cd_total,Cl_total,Aref,Cd_sent,Cd_scham] = single_plate_quasi_cd(phi, theta, S, Ai, m_r, Ri, Tw, Vi, A, Ar, Alpha, f, flag_axis, v_sbf)
% provides the capability of computing drag-coefficient for panels with
% different angles of rotation in the body frame 
if strcmp(flag_axis,'z') == 1
    u = [cos(phi).*cos(theta); cos(phi).*sin(theta); sin(phi)];
elseif strcmp(flag_axis,'y') == 1
    u = [cos(phi).*cos(theta); sin(phi); cos(phi).*sin(theta)];   %*ones(1,numel(theta))
elseif strcmp(flag_axis,'x') == 1
    u = [sin(phi); cos(phi).*cos(theta); cos(phi).*sin(theta);];
end
if ~isempty(v_sbf)
    u = v_sbf;
end
Cd_scham = 0;  %%%%%%%%%%%%%% 0
Aref = 0;
for j=1:numel(Ai(1,:))
    Y = u(1)*Ai(1,j) + u(2)*Ai(2,j) + u(3)*Ai(3,j);
%     P = exp(-Y.*Y.*(S*S))./S;
    Alpha_scham = 3.6*m_r(j)*Y/(1+m_r(j))^2;
    r_sch = sqrt(1+Alpha_scham*(3*Ri*Tw/Vi^2 -1));     % Koppenwallner's velocity ratio
    Y(Y<0) = 0;
    Aref = Aref + A(j)*Y;
    Cd_scham = Cd_scham + 2*A(j)*Y.*(1+r_sch.*(2*Y.^2-1)); %+ A(j)*P/sqrt(pi);
end
Cd_scham = Cd_scham/Ar;
for j = 1:numel(Ai(1,:))
    Y(j) = u(1)*Ai(1,j) + u(2)*Ai(2,j) + u(3)*Ai(3,j);
    ul = cross(cross(u,Ai(:,j)), u);
    ul = ul/norm(ul);
    l = ul(1)*Ai(1,j) + ul(2)*Ai(2,j) + ul(3)*Ai(3,j);
    [Cd(j,:),Cl_ads(j,:)] = sentman(Vi, A(j), Ar,Y(j),Ri,Alpha,Tw,S,l);
end
Cd_sent = sum(Cd)+ 0;  %%%%%%%%%%%%%%%%%%%%%%
Cd_total = f*Cd_sent + (1-f)*Cd_scham;
Cl_total = f*Cl_ads;