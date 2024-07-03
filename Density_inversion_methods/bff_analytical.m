%% calculate the BFF coefficients
function [An_total, Bn_total] = bff_analytical(A, Ar, Ai, S, r_ads, r_s, frac, order_body,flag_axis,flag_method,flag_diff)

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
        An_ads = fourier_plate_dria(A, Ar, S, r_ads,Ci, order_body);
%         r_s = r_s*ones(1,numel(A));
        An_s = fourier_plate_dria(A, Ar, S, r_s, Ci,order_body);
        An_plates = frac*An_ads + (1-frac)*An_s;
    case 'diff'
        switch flag_diff
            case 's'
                r_ads = r_ads*ones(1,numel(A));
%                 r_s = r_s*ones(1,numel(A));
%                 An_ads = fourier_plate_dria_diff(A, Ar, S, r_ads,Ci, order_body,flag_diff);
%                 An_s = fourier_plate_dria_diff(A, Ar, S, r_s, Ci,order_body,flag_diff);
                An_ads1 = fourier_plate_dria(A, Ar, S-0.05, r_ads,Ci, order_body);
                An_ads2 = fourier_plate_dria(A, Ar, S+0.05, r_ads,Ci, order_body);
                An_ads = (An_ads2-An_ads1)/0.1;
                An_s1 = fourier_plate_dria(A, Ar, S-0.05, r_s,Ci, order_body);
                An_s2 = fourier_plate_dria(A, Ar, S+0.05, r_s,Ci, order_body);
                An_s = (An_s2-An_s1)/0.1;                
                An_plates = frac*An_ads + (1-frac)*An_s;
            case 'r_ads'
                r_ads = r_ads*ones(1,numel(A));
                An_ads1 = fourier_plate_dria(A, Ar, S, r_ads-0.05,Ci, order_body);
                An_ads2 = fourier_plate_dria(A, Ar, S, r_ads+0.05,Ci, order_body);
                An_plates = frac*(An_ads2-An_ads1)/0.1;
            case 'r_s'
                r_s = r_s*ones(1,numel(A));
                An_s1 = fourier_plate_dria(A, Ar, S, r_s-0.05,Ci, order_body);
                An_s2 = fourier_plate_dria(A, Ar, S, r_s+0.05,Ci, order_body);
                An_plates = (1-frac)*(An_s2-An_s1)/0.1;
            case 'f'
                r_ads = r_ads*ones(1,numel(A));
                An_ads = fourier_plate_dria(A, Ar, S, r_ads,Ci, order_body);
%                 r_s = r_s*ones(1,numel(A));
                An_s = fourier_plate_dria(A, Ar, S, r_s, Ci,order_body);
                An_plates = An_ads-An_s;
        end
end



order_vec_body = [0:order_body]';

cos_delta = cos(order_vec_body.*delta);
sin_delta = sin(order_vec_body.*delta);

An_total = sum(An_plates.*cos_delta,2);
Bn_total = sum(An_plates.*sin_delta,2);