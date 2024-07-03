function [area_vec, Area_plates, spec_ref, diff_ref, spec_ref_ir, diff_ref_ir, M_s, mass] = spire_properties_simple(version_spire, area)

switch version_spire
    case '3.0'
        area_vec(:,1) = [0;0;1]; % white aeroglaze
        area_vec(:,2) = [0;0;-1]; % PCB
        area_vec(:,3) = [0;1;0]; % Anodized Aluminum
        area_vec(:,4) = [0;-1;0]; % white aeroglaze
        area_vec(:,5) = [1;0;0]; % PCB
        area_vec(:,6) = [-1;0;0]; % Anodized Aluminum

        area_vec(:,7) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half (PCB)
        area_vec(:,8) = sqrt(2)/2*[ -1;1; 0];               % front solar panel half (Solar cell)

        Area_plates = [0.1*0.1, 0.1*0.1, 0.338*0.1, 0.338*0.1, 0.338*0.1, 0.338*0.1, 0.32*0.15*2, 0.32*0.15*2]/area;
        
        
        spec_ref = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5];
        diff_ref = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2];
        spec_ref_ir = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5];
        diff_ref_ir = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2];
        M_s = [27, 27, 27, 27, 27, 27, 100, 100, 100, 100];
        mass = 4.455;
    case '3.3'
        area_vec(:,1) = [0;0;1]; % white aeroglaze
        area_vec(:,2) = [0;0;-1]; % PCB
        area_vec(:,3) = [0;1;0]; % Anodized Aluminum
        area_vec(:,4) = [0;-1;0]; % white aeroglaze
        area_vec(:,5) = [1;0;0]; % PCB
        area_vec(:,6) = [-1;0;0]; % Anodized Aluminum

        area_vec(:,7) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half (PCB)
        area_vec(:,8) = sqrt(2)/2*[ -1;1; 0];               % front solar panel half (Solar cell)

        Area_plates = [0.1*0.1, 0.1*0.1, 0.338*0.1, 0.338*0.1, 0.338*0.1, 0.338*0.1, 0.32*0.15*2, 0.32*0.15*2]/area;
        
        
        spec_ref = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5];
        diff_ref = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2];
        spec_ref_ir = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5];
        diff_ref_ir = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2];
        M_s = [27, 27, 27, 27, 27, 27, 100, 100, 100, 100];
        mass  = 4.933;
    case '3.7'
        area_vec(:,1) = [0;0;1]; % white aeroglaze
        area_vec(:,2) = [0;0;-1]; % PCB
        area_vec(:,3) = [0;1;0]; % Anodized Aluminum
        area_vec(:,4) = [0;-1;0]; % white aeroglaze
        area_vec(:,5) = [1;0;0]; % PCB
        area_vec(:,6) = [-1;0;0]; % Anodized Aluminum

        area_vec(:,7) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half (PCB)
        area_vec(:,8) = sqrt(2)/2*[ -1;1; 0];               % front solar panel half (Solar cell)

        Area_plates = [0.1*0.1, 0.1*0.1, 0.338*0.1, 0.338*0.1, 0.338*0.1, 0.338*0.1, 0.32*0.224*2, 0.32*0.224*2]/area;
        
        
        spec_ref = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5];
        diff_ref = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2];
        spec_ref_ir = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5];
        diff_ref_ir = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2];
        M_s = [27, 27, 27, 27, 27, 27, 100, 100, 100, 100];
        mass = 5.1;   
    case '4.11'
        area_vec(:,1) = [0;0;1]; % white aeroglaze
        area_vec(:,2) = [0;0;-1]; % PCB
        area_vec(:,3) = [0;1;0]; % Anodized Aluminum
        area_vec(:,4) = [0;-1;0]; % white aeroglaze
        area_vec(:,5) = [1;0;0]; % PCB
        area_vec(:,6) = [-1;0;0]; % Anodized Aluminum

        area_vec(:,7) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half (PCB)
        area_vec(:,8) = sqrt(2)/2*[ -1;1; 0];               % front solar panel half (Solar cell)

        Area_plates = [0.1*0.1, 0.1*0.1, 0.338*0.1, 0.338*0.1, 0.338*0.1, 0.338*0.1, 0.32*0.224*2, 0.32*0.224*2]/area;
        
        
        spec_ref = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5];
        diff_ref = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2];
        spec_ref_ir = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5];
        diff_ref_ir = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2];
        M_s = [27, 27, 27, 27, 27, 27, 100, 100, 100, 100];
        mass = 5.1;        
end