function hasdm_mat = hasdm_initialize(n_days, doy, hasdm_models)
hasdm_mat = [];
for ii = 1:numel(hasdm_models)
    ind_day = [];
    load(hasdm_models{ii}, 'den_mat')
    for nn = 1:n_days
        doy_hasdm = doy+nn-1;
        ind_day  = [ind_day;find(den_mat(:,2) == doy_hasdm)];
    end
    hasdm_mat = [hasdm_mat;den_mat(ind_day,:)];
    clear den_mat 
end