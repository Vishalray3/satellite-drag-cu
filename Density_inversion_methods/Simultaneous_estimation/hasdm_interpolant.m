function F_hasdm = hasdm_interpolant(n_days, doy, den_mat)


ind_day = [];
for nn = 1:n_days
    doy_hasdm = doy+nn-1;
    ind_day  = [ind_day;find(den_mat(:,2) == doy_hasdm)];
end
hasdm_mat = den_mat(ind_day,:);


vec_alt  = unique(hasdm_mat(:,7));
[vec_njd, ind_njd_all]  = unique(hasdm_mat(:,6));
vec_long = unique(hasdm_mat(:,10));
vec_lat  = unique(hasdm_mat(:,9));
hasdm_mat(:,11) = log(hasdm_mat(:,11));
ind_long = [1:numel(vec_long)];
den_grid = zeros(numel(vec_alt), numel(vec_njd),numel(vec_long), numel(vec_lat));

N_njd = numel(vec_njd);
N_lat = numel(vec_lat);
N_alt = numel(vec_alt);
for ii_alt = 1:N_alt
    ind_alt = find(hasdm_mat(:,7) == vec_alt(ii_alt));
    for ii_njd = 1:N_njd
        ind_njd = find(hasdm_mat(:,6) == vec_njd(ii_njd));
        for ii_lat = 1:N_lat
            ind_alt_njd = intersect(ind_alt,ind_njd);
            hasdm_small = hasdm_mat(ind_alt_njd,:);
            
            ind_lat = hasdm_small(:,9) == vec_lat(ii_lat);
            hasdm_small = hasdm_small(ind_lat,:);
            
            den_ind = interp1(hasdm_small(:,10), hasdm_small(:,11), vec_long,'linear', 'extrap');
            den_grid(ii_alt,ii_njd,ind_long,ii_lat) = den_ind;
        end
    end
end
jd_ref = GREGORIANtoJD_vector(1950,1,0,0,0,0);

% njd calculation to be uniform between SET HASDM and Eric's data
yearh = hasdm_mat(ind_njd_all,1);
dayofyh = hasdm_mat(ind_njd_all,2);
hrh = hasdm_mat(ind_njd_all,3);
minh = hasdm_mat(ind_njd_all,4);
sech = hasdm_mat(ind_njd_all,5);
jdutc = GREGORIANtoJD_vector(yearh, dayofyh, hrh, minh, sech) ;
vec_njd_correct = jdutc - jd_ref;

gridVecs = {vec_alt,vec_njd_correct,vec_long,vec_lat};
F_hasdm = griddedInterpolant(gridVecs,den_grid);