%% function
function [P_total, Q_total, An_total, Bn_total] = torque_gsim_inversion(parameters, flag_method,flag_diff)
Ar = parameters.Ar;
A = parameters.Area_plates;
order_body = parameters.Order_b;
Ai = parameters.area_vec;
flag_axis = parameters.flag_axis;
Cp = parameters.Cp;

S = parameters.S;


frac =  parameters.frac; %Kl.*P_o./(1+Kl.*P_o);%
r_ads = parameters.r_ads;
r_s = parameters.r_s;

if strcmp(flag_axis,'y')
    Ci = sqrt(Ai(1,:).^2+Ai(3,:).^2);
    delta = atan2(Ai(3,:), Ai(1,:));
elseif strcmp(flag_axis, 'z')
    Ci = sqrt(Ai(1,:).^2+Ai(2,:).^2);
    delta = atan2(Ai(1,:), Ai(2,:));
end

switch flag_method
    case 'coeff'
        r_ads = r_ads*ones(1,numel(A));
        An_ads = fourier_plate_dria(A, Ar, S, r_ads,Ci, order_body+1);
%         r_s = r_s*ones(1,numel(A));
        An_s = fourier_plate_dria(A, Ar, S, r_s, Ci,order_body+1);
        An_plates = frac*An_ads + (1-frac)*An_s;
    case 'diff'
        switch flag_diff
            case 's'
                r_ads = r_ads*ones(1,numel(A));
%                 r_s = r_s*ones(1,numel(A));
%                 An_ads = fourier_plate_dria_diff(A, Ar, S, r_ads,Ci, order_body,flag_diff);
%                 An_s = fourier_plate_dria_diff(A, Ar, S, r_s, Ci,order_body,flag_diff);
                An_ads1 = fourier_plate_dria(A, Ar, S-0.05, r_ads,Ci, order_body+1);
                An_ads2 = fourier_plate_dria(A, Ar, S+0.05, r_ads,Ci, order_body+1);
                An_ads = (An_ads2-An_ads1)/0.1;
                An_s1 = fourier_plate_dria(A, Ar, S-0.05, r_s,Ci, order_body+1);
                An_s2 = fourier_plate_dria(A, Ar, S+0.05, r_s,Ci, order_body+1);
                An_s = (An_s2-An_s1)/0.1;                
                An_plates = frac*An_ads + (1-frac)*An_s;
            case 'r_ads'
                r_ads = r_ads*ones(1,numel(A));
                An_ads1 = fourier_plate_dria(A, Ar, S, r_ads-0.05,Ci, order_body+1);
                An_ads2 = fourier_plate_dria(A, Ar, S, r_ads+0.05,Ci, order_body+1);
                An_plates = frac*(An_ads2-An_ads1)/0.1;
            case 'r_s'
                r_s = r_s*ones(1,numel(A));
                An_s1 = fourier_plate_dria(A, Ar, S, r_s-0.05,Ci, order_body+1);
                An_s2 = fourier_plate_dria(A, Ar, S, r_s+0.05,Ci, order_body+1);
                An_plates = (1-frac)*(An_s2-An_s1)/0.1;
            case 'f'
                r_ads = r_ads*ones(1,numel(A));
                An_ads = fourier_plate_dria(A, Ar, S, r_ads,Ci, order_body+1);
%                 r_s = r_s*ones(1,numel(A));
                An_s = fourier_plate_dria(A, Ar, S, r_s, Ci,order_body+1);
                An_plates = An_ads-An_s;
        end
end


%% Torque Fourier calculation
order_vec_body = [0:order_body+1]';

cos_delta = cos(order_vec_body.*delta);
sin_delta = sin(order_vec_body.*delta);
An_total = An_plates.*cos_delta;
Bn_total = An_plates.*sin_delta;

if strcmp(flag_axis, 'z')
    F_vec = [An_total(2,:); Bn_total(2,:); zeros(1, numel(delta))];
    P_total(:,1)  = sum(cross(Cp, F_vec, 1),2);
    Q_total(:,1) = [0;0;0];
    
    if order_body > 0
        F_vec = [(2*An_total(1,:) + An_total(3,:))/2; Bn_total(3,:)/2; zeros(1, numel(delta))];
        P_total(:,2)  = sum(cross(Cp, F_vec, 1),2);
        F_vec = [Bn_total(3,:)/2; (2*An_total(1,:) - An_total(3,:))/2; zeros(1, numel(delta))];
        Q_total(:,2)  = sum(cross(Cp, F_vec, 1),2);
        if order_body > 1
            for n = 3:order_body+1
                F_vec = [(An_total(n+1,:) + An_total(n-1,:))/2; (Bn_total(n+1,:) - Bn_total(n-1,:))/2; zeros(1, numel(delta))];
                P_total(:,n)  = sum(cross(Cp, F_vec, 1),2);
                F_vec = [(Bn_total(n+1,:) + Bn_total(n-1,:))/2; (An_total(n-1,:) - An_total(n+1,:))/2; zeros(1, numel(delta))];
                Q_total(:,n)  = sum(cross(Cp, F_vec, 1),2);
            end
            
        end
        
    end
    
end

end