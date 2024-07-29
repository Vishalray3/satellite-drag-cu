function hasdm_mat = hasdm_initialize(n_days, doy, den_mat_list)
hasdm_mat = [];
for ii = 1:numel(den_mat_list)
    ind_day = [];
    den_mat = den_mat_list{ii};
    for nn = 1:n_days
        doy_hasdm = doy+nn-1;
        ind_day  = [ind_day;find(den_mat(:,2) == doy_hasdm)];
    end
    hasdm_mat = [hasdm_mat;den_mat(ind_day,:)];
end