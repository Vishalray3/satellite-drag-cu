%% read fes ocean tides data
clc
clearvars
ndeg = 10;
nord = 10;
fes = readmatrix('fes2004_Cnm-Snm.txt');

doo = unique(fes(:,1));

Cnmp = fes(:,5);
Snmp = fes(:,6);
Cnmm = fes(:,7);
Snmm = fes(:,8);

for jj = 2:ndeg
    vec_deg = find(fes(:,3) == jj);
    for kk = 0:nord
        lin_ind = jj*(jj+1)/2 + kk - 2;
        vec_ord = find(fes(:,4) == kk);
        vec_do = intersect(vec_deg,vec_ord);
        doo_do = fes(vec_do,1);
        Cnmp_do = Cnmp(vec_do);
        Snmp_do = Snmp(vec_do);
        Cnmm_do = Cnmm(vec_do);
        Snmm_do = Snmm(vec_do);        
        deg_do(lin_ind) = jj;
        ord_do(lin_ind) = kk;
        for ii = 1:numel(doo)
           ind_do = find(doo_do == doo(ii));
           if ~isempty(ind_do)
               Cnmp_ot(ii, lin_ind)= Cnmp_do(ind_do);
               Snmp_ot(ii, lin_ind)= Snmp_do(ind_do);
               Cnmm_ot(ii, lin_ind)= Cnmm_do(ind_do);
               Snmm_ot(ii, lin_ind)= Snmm_do(ind_do);               
           end
        end
    end
end
doo = doo*1000;
doo_str = num2str(doo);
doo_coeff(:,1) = [zeros(8,1); ones(4,1); 2*ones(5,1); 4];
for ii = 2:numel(doo_str(1,:))
doo_coeff(:,ii) = str2num(doo_str(:,ii))-5;
end

save('fes2010_oceantide','Cnmp_ot','Snmp_ot','Cnmm_ot','Snmm_ot','doo_coeff','deg_do','ord_do')