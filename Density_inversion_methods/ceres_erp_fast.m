classdef ceres_erp_fast
    % Written by Vishal Ray 03/2022
    % Ref: Doornbos et al. air density models derived from multiple
    % satellite observations
    properties
        sunvector_j2000        % Sun coordinates in J2000 frame (m)
        nrm_eci                % panel normal vectors in eci frame
        plt_area               % panel areas (m^2)
        scmass                 % satellite mass (kg)
        refl_vis_spec          % visible specular reflectivity coefficient
        refl_vis_diff          % visible diffusive reflectivity coefficient
        refl_ir_spec           % infrared specular reflectivity coefficient
        refl_ir_diff           % infrared diffusive reflectivity coefficient
        ECEFtoJ2000            % rotation matrix from ECEF to J2000 frame
        Psw                    % outgoing shortwave radiation pressure (flux/c)
        Plw                    % outgoing longwave radiation pressure (flux/c)
        lam                    % longitude of grid (degrees)
        phi                    % latitude of grid (degrees)
        Cr_erp                 % scaling factor
        Re                     % radius of Earth
    end
    
    methods
        function obj =  ceres_erp_fast(sunvector_j2000,nrm_eci,plt_area,scmass,refl_vis_spec,refl_vis_diff,refl_ir_spec,refl_ir_diff,ECEFtoJ2000, ...
                Psw, Plw, lam, phi, Cr_erp, Re)
            obj.sunvector_j2000     = sunvector_j2000;
            obj.nrm_eci    = nrm_eci;
            obj.plt_area   = plt_area;
            obj.scmass     = scmass;
            obj.refl_vis_spec  = refl_vis_spec;
            obj.refl_vis_diff  = refl_vis_diff;
            obj.refl_ir_spec   = refl_ir_spec;
            obj.refl_ir_diff   = refl_ir_diff;
            obj.ECEFtoJ2000    = ECEFtoJ2000;
            obj.Psw = Psw;
            obj.Plw = Plw;
            obj.lam = lam;
            obj.phi = phi;
            obj.Cr_erp = Cr_erp;
            obj.Re = Re;
        end
         function [a, F_ar, F_av, F_aPar] = compAccelPar(obj, X_state,~)
            
            
            phi_rad = obj.phi*pi/180;
            lam_rad = obj.lam*pi/180;
            [a, F_ar] = ceres_erp_mex(obj.sunvector_j2000', obj.nrm_eci, obj.plt_area, obj.scmass, obj.refl_vis_spec, obj.refl_vis_diff, obj.refl_ir_spec, ...
                obj.refl_ir_diff, obj.ECEFtoJ2000, obj.Psw, obj.Plw, lam_rad, phi_rad, obj.Cr_erp, obj.Re, X_state(1:3)', numel(obj.phi), numel(obj.lam), ...
                numel(obj.plt_area));

            F_av = zeros(3,3);
            F_aPar = a/obj.Cr_erp;  
            
         end
        
        function a = compAccel(obj, X_state,~)
            
            phi_rad = obj.phi*pi/180;
            lam_rad = obj.lam*pi/180;
            a = ceres_erp_mex(obj.sunvector_j2000', obj.nrm_eci, obj.plt_area, obj.scmass, obj.refl_vis_spec, obj.refl_vis_diff, obj.refl_ir_spec, ...
                obj.refl_ir_diff, obj.ECEFtoJ2000, obj.Psw, obj.Plw, lam_rad, phi_rad, obj.Cr_erp, obj.Re, X_state(1:3)', numel(obj.phi), numel(obj.lam), ...
                numel(obj.plt_area));
           
        end
    end
end