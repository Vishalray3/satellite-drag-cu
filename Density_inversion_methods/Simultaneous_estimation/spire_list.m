%% Lemur versions

function version_spire = spire_list(sat_ID)

switch sat_ID
    case {'FM046','FM047', 'FM048','FM049','FM050','FM051','FM052','FM053','FM054','FM055','FM056','FM057','FM059','FM060','FM061','FM062','FM063','FM064',...
            'FM065','FM066','FM067','FM075','FM078','FM079','FM080'}
        version_spire = '3.0';
        
    case {'FM081','FM082','FM083','FM084','FM085','FM086','FM087','FM088','FM089','FM090'}
        version_spire = '3.3';
        
    case {'FM091','FM092','FM093','FM094','FM095','FM096','FM097','FM098','FM099','FM100','FM101','FM102','FM103','FM104','FM105','FM106','FM107','FM108',...
'FM109', 'FM110','FM113','FM115','FM117','FM118','FM119','FM120','FM122','FM124','FM125','FM126'}
        version_spire = '3.7';
        
    case {'FM128','FM129','FM132','FM133','FM134','FM135','FM143','FM148','FM150'}
        version_spire = '4.11';    
end