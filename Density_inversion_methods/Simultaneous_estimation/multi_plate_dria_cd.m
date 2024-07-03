function [Cd_total,Aref] = multi_plate_dria_cd(phi, theta, S, Ai, m_r, Ri, Tw, Vi, A, Ar, Alpha, f, flag_axis, v_sbf)
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
Cd_s = 0;  %%%%%%%%%%%%%% 0
for j = 1:numel(Ai(1,:))
    Y = u(1,j)*Ai(1,j) + u(2,j)*Ai(2,j) + u(3,j)*Ai(3,j);
    Alpha_s = 3.6*m_r(j)/(1+m_r(j))^2;  % 3.6*m_r(j)*Y/(1+m_r(j))^2
    Cd_s(j,:) = sentman(Vi, A(j), Ar,Y,Ri,Alpha_s,Tw,S);
end

Cd_s = sum(Cd_s);
for j = 1:numel(Ai(1,:))
    Y(j) = u(1,j)*Ai(1,j) + u(2,j)*Ai(2,j) + u(3,j)*Ai(3,j);
    Cd(j,:) = sentman(Vi, A(j), Ar,Y(j),Ri,Alpha,Tw,S);
end
Aref = sum(A(Y>0).*Y(Y>0));
Cd_ads = sum(Cd)+ 0;  %%%%%%%%%%%%%%%%%%%%%%
Cd_total = f*Cd_ads + (1-f)*Cd_s;