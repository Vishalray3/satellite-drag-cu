function den_grid = hasdm_den_grid(hasdm_mat)
%#codegen
vec_alt  = unique(hasdm_mat(:,7));
vec_njd  = unique(hasdm_mat(:,6));
vec_long = unique(hasdm_mat(:,10));
vec_lat  = unique(hasdm_mat(:,9));
%%%%%%%%%%
%         if numel(vec_long) > 24
%             vec_long = (vec_long(1:2:end) + vec_long(2:2:end))/2;   %%% Verify for each date
%         end
%%%%%%%%%
den_grid = zeros(numel(vec_alt), numel(vec_njd),numel(vec_long), numel(vec_lat));
for ii_alt = 1:numel(vec_alt)
    ind_alt = find(hasdm_mat(:,7) == vec_alt(ii_alt));
    for ii_njd = 1:numel(vec_njd)
        ind_njd = find(hasdm_mat(:,6) == vec_njd(ii_njd));
        for ii_long = 1:numel(vec_long)
            ind_long = find(abs(hasdm_mat(:,10) - vec_long(ii_long)) < 1);
            %                     ind_long = find(hasdm_mat(:,10) == vec_long(ii_long));
            for ii_lat = 1:numel(vec_lat)
                ind_lat = find(hasdm_mat(:,9) == vec_lat(ii_lat));
                ind_den = intersect(intersect(intersect(ind_alt,ind_njd), ind_long), ind_lat);
                den_grid(ii_alt,ii_njd,ii_long,ii_lat) = hasdm_mat(ind_den(1),11);
            end
        end
    end
end