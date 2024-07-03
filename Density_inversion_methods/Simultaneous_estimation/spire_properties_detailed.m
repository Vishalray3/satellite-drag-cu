function [area_vec, Area_plates, spec_ref, diff_ref, spec_ref_ir, diff_ref_ir, M_s, mass] = spire_properties_detailed(version_spire)

switch version_spire
    case '3.0'
        area_vec(:,1) = [1;0;0];
        area_vec(:,2) = [-1;0;0];
        area_vec(:,3) = [0;1;0];
        area_vec(:,4) = [0;-1;0];
        area_vec(:,5) = [0;0;1];
        area_vec(:,6) = [0;0;-1];
        area_vec(:,7) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half
        area_vec(:,8) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half
        area_vec(:,9) = sqrt(2)/2*[-1; 1; 0];              % Aft solar panel
        
        % Area_plates = [0.1053314*0.3388226, 0.1053314*0.3388226, 0.105186*0.3388226 + 3*0.005334*0.074549, 0.105186*0.3388226 + 3*0.005334*0.074549, ...
        %     0.105186*0.1053314, 0.105186*0.1053314, 0.150*0.3222, 0.150*0.3222,
        %     2*0.150*0.3222]/area; eric's code
        Area_plates = [0.1053314*0.3388226, 0.1053314*0.3388226, 0.105186*0.3388226 + 3*0.005334*0.074549, 0.105186*0.3388226 + 3*0.005334*0.074549, ...
            0.105186*0.1053314 + 0.075*0.105, 0.105186*0.1053314 + 0.075*0.105, 0.150*0.3222, 0.150*0.3222, 2*0.150*0.3222]/area;
        
        spec_ref = [0.2, 0.05, 0.05, 0.2, 0.2, 0.2, 0.2, 0.05, 0.05];
        diff_ref = [0.4, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3];
        spec_ref_ir = [0.2, 0.05, 0.05, 0.2, 0.2, 0.2, 0.2, 0.05, 0.05];
        diff_ref_ir = [0.4, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3];
        M_s = [26.981538, 60.0843, 60.0843, 26.981538, 26.981538, 26.981538, 26.981538, 60.0843, 60.0843];
        mass = 4.455;
    case '3.3'
        area_vec(:,1) = [1;0;0]; % white aeroglaze
        area_vec(:,2) = [1;0;0]; % PCB
        area_vec(:,3) = [1;0;0]; % Anodized Aluminum
        area_vec(:,4) = [-1;0;0]; % solar cell
        area_vec(:,5) = [-1;0;0]; % PCB
        
        area_vec(:,6) = [0;1;0];  % solar cell
        area_vec(:,7) = [0;1;0]; % PCB 
        area_vec(:,8) = [0;-1;0];  % white aeroglaze
        area_vec(:,9) = [0;-1;0]; % Ano Alum 
        
        area_vec(:,10) = [0;0;1]; % pcb
        area_vec(:,11) = [0;0;-1]; % pcb
        
        area_vec(:,12) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half (PCB)
        area_vec(:,13) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half (solar cell)
        area_vec(:,14) = sqrt(2)/2*[-1; 1; 0];               % Aft solar panel (pcb)
        area_vec(:,15) = sqrt(2)/2*[-1; 1; 0];               % Aft solar panel (solar cell)
        
        
        Area_plates = [0.27*0.105, 0.032*0.105, 0.032*0.105, 0.069*0.0385*7, 0.0165, ....
            0.069*0.0385*7, 0.0165, 0.27*0.105 + 0.007*0.7*3, 0.064*0.105,...
            0.105*0.105,0.105*0.105, ...
            0.0165*2 + 0.0750*0.3222*2, 0.069*0.0385*7*2, 0.0165*4, 0.069*0.0385*7*4]/area;
        % the knob like projection - need to use a cylinder to model this
        spec_ref = [0.3, 0.12, 0.5, 0.07, 0.12, 0.07, 0.12, 0.3, 0.5, 0.12, 0.12, 0.12, 0.07, 0.12, 0.07]; 
        diff_ref = [0.46, 0.05, 0.2, 0.02, 0.05, 0.02, 0.05, 0.46, 0.2, 0.05, 0.05, 0.05, 0.02, 0.05, 0.02];
        spec_ref_ir = [0.3, 0.12, 0.5, 0.07, 0.12, 0.07, 0.12, 0.3, 0.5, 0.12, 0.12, 0.12, 0.07, 0.12, 0.07]; 
        diff_ref_ir = [0.46, 0.05, 0.2, 0.02, 0.05, 0.02, 0.05, 0.46, 0.2, 0.05, 0.05, 0.05, 0.02, 0.05, 0.02];
        M_s = [88.1, 588, 26.9, 60.1, 588, 60.1, 588, 88.1, 26.9, 88.1, 88.1, 88.1, 60.1, 88.1, 60.1];
        mass  = 4.933;
    case '3.7'
        area_vec(:,1) = [1;0;0]; % white aeroglaze
        area_vec(:,2) = [1;0;0]; % PCB
        area_vec(:,3) = [1;0;0]; % Anodized Aluminum
        area_vec(:,4) = [-1;0;0]; % white aeroglaze
        area_vec(:,5) = [-1;0;0]; % PCB
        area_vec(:,6) = [-1;0;0]; % Anodized Aluminum
        
        area_vec(:,7) = [0;1;0];  % white aeroglaze
        area_vec(:,8) = [0;1;0]; % Ano Alum 
        area_vec(:,9) = [0;-1;0];  % white aeroglaze
        area_vec(:,10) = [0;-1;0]; % Ano Alum 
        
        area_vec(:,11) = [0;0;1]; % pcb
        area_vec(:,12) = [0;0;-1]; % pcb
        
        area_vec(:,13) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half (PCB)
        area_vec(:,14) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half (Solar cell)
        area_vec(:,15) = sqrt(2)/2*[-1; 1; 0];              % Aft solar panel (PCB)
        area_vec(:,16) = sqrt(2)/2*[-1; 1; 0];              % Aft solar panel (Solar cell)
        % Area_plates = [0.1053314*0.3388226, 0.1053314*0.3388226, 0.105186*0.3388226 + 3*0.005334*0.074549, 0.105186*0.3388226 + 3*0.005334*0.074549, ...
        %     0.105186*0.1053314, 0.105186*0.1053314, 0.150*0.3222, 0.150*0.3222,
        %     2*0.150*0.3222]/area; eric's code
        
        Area_plates = [0.27*0.105, 0.032*0.105, 0.032*0.105, 0.27*0.105, 0.032*0.105, 0.032*0.105,...
             0.27*0.105 + 0.007*0.7*3, 0.27*0.105 + 0.007*0.7*3, 0.27*0.105 + 0.007*0.7*3, 0.27*0.105 + 0.007*0.7*3,  ...
             0.105*0.105, 0.105*0.105,...
            0.0165*2 + 0.0750*0.3222*4, 0.069*0.0385*7*2 ,  0.0165*6, 0.069*0.0385*7*6]/area;
        
        
        % the knob like projection - need to use a cylinder to model this
        spec_ref = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5, 0.3, 0.5, ...
                    0.12, 0.12, ...
                    0.12, 0.07, 0.12, 0.07];
        diff_ref = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2,0.46, 0.2,...
                    0.05, 0.05,...
                    0.05, 0.02, 0.05, 0.02];
        spec_ref_ir = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5, 0.3, 0.5, ...
                    0.12, 0.12, ...
                    0.12, 0.07, 0.12, 0.07];
        diff_ref_ir = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2,0.46, 0.2,...
                    0.05, 0.05,...
                    0.05, 0.02, 0.05, 0.02];
        M_s = [88.1, 588, 26.9, 88.1, 588, 26.9, 88,1, 26.9, 88.1, 26.9, 588, 588, 588, 60.1, 588, 60.1];
        mass = 5.1;
    case '4.11'
        area_vec(:,1) = [1;0;0]; % white aeroglaze
        area_vec(:,2) = [1;0;0]; % PCB
        area_vec(:,3) = [1;0;0]; % Anodized Aluminum
        area_vec(:,4) = [-1;0;0]; % white aeroglaze
        area_vec(:,5) = [-1;0;0]; % PCB
        area_vec(:,6) = [-1;0;0]; % Anodized Aluminum
        
        area_vec(:,7) = [0;1;0];  % white aeroglaze
        area_vec(:,8) = [0;1;0]; % Ano Alum 
        area_vec(:,9) = [0;-1;0];  % white aeroglaze
        area_vec(:,10) = [0;-1;0]; % Ano Alum 
        
        area_vec(:,11) = [0;0;1]; % pcb
        area_vec(:,12) = [0;0;-1]; % pcb
        
        area_vec(:,13) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half (PCB)
        area_vec(:,14) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half (Solar cell)
        area_vec(:,15) = sqrt(2)/2*[-1; 1; 0];              % Aft solar panel (PCB)
        area_vec(:,16) = sqrt(2)/2*[-1; 1; 0];              % Aft solar panel (Solar cell)
        % Area_plates = [0.1053314*0.3388226, 0.1053314*0.3388226, 0.105186*0.3388226 + 3*0.005334*0.074549, 0.105186*0.3388226 + 3*0.005334*0.074549, ...
        %     0.105186*0.1053314, 0.105186*0.1053314, 0.150*0.3222, 0.150*0.3222,
        %     2*0.150*0.3222]/area; eric's code
        
        Area_plates = [0.27*0.105, 0.032*0.105, 0.032*0.105 + 0.039*0.025, 0.27*0.105, 0.032*0.105, 0.032*0.105 + 0.039*0.025,...
             0.27*0.105 + 0.007*0.7*3, 0.27*0.105 + 0.007*0.7*3 + 0.039*0.025, 0.27*0.105 + 0.007*0.7*3, 0.27*0.105 + 0.007*0.7*3 + 0.039*0.025,  ...
             0.105*0.105, 0.105*0.105,...
            0.0165*2 + 0.0750*0.3222*4, 0.069*0.0385*7*2 ,  0.0165*6, 0.069*0.0385*7*6]/area;
        
        
        % the knob like projection - need to use a cylinder to model this
        spec_ref = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5, 0.3, 0.5, ...
                    0.12, 0.12, ...
                    0.12, 0.07, 0.12, 0.07];
        diff_ref = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2,0.46, 0.2,...
                    0.05, 0.05,...
                    0.05, 0.02, 0.05, 0.02];
        spec_ref_ir = [0.3, 0.12, 0.5, 0.3, 0.12, 0.5, ...
                    0.3, 0.5, 0.3, 0.5, ...
                    0.12, 0.12, ...
                    0.12, 0.07, 0.12, 0.07];
        diff_ref_ir = [0.46, 0.05, 0.2, 0.46, 0.05, 0.2,...
                    0.46, 0.2,0.46, 0.2,...
                    0.05, 0.05,...
                    0.05, 0.02, 0.05, 0.02];
        M_s = [88.1, 588, 26.9, 88.1, 588, 26.9, 88,1, 26.9, 88.1, 26.9, 588, 588, 588, 60.1, 588, 60.1];
        mass = 5.2;        
end