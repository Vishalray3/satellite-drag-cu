%% function
function a_srp = srp_nadir_fourier(angle_step,parameters,sdel, order_fourier,trig)

slam = angle_step;
n_hat = parameters.area_vec;
area = parameters.Area_plates;
rho_spec = parameters.rho_spec;
rho_diff = parameters.rho_diff;

sun_vec = [cos(sdel)*sin(slam);sin(sdel);cos(sdel)*cos(slam)]; %X_rel/r_rel;
sun_vec_mat = repmat(sun_vec,1,numel(area));
cos_theta = dot(sun_vec_mat, n_hat);
%             sp_angle = acosd(cos_theta);
cos_theta(cos_theta<0) = 0;
e_coeff = sum(area.*cos_theta.*(1-rho_spec));
n_coeff_mat = repmat(area.*rho_spec.*cos_theta.^2 + area.*rho_diff.*cos_theta/3,3,1);
n_comp = sum(n_coeff_mat.*n_hat,2);
a_srp = (e_coeff*sun_vec + 2*n_comp);





order_vec_orbit = [0:order_fourier];
if trig ==1
    a_srp = repmat(a_srp, 1, order_fourier+1).*repmat(cos(order_vec_orbit*slam),3,1);
elseif trig ==2
    a_srp = repmat(a_srp, 1, order_fourier+1).*repmat(sin(order_vec_orbit*slam),3,1); 
end
end