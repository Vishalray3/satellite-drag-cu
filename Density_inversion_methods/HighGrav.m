classdef HighGrav
    properties
        et                                              % ephemeris time
        deg                                             % degree of model
        ord                                             % order of model
        mu_e                                            % grav coeff in m3/s2  
        Re                                              % radius of earth in m 
        Cbar                                            % Normalized grav. coeff. in 2d matrix, see GravCoeffArrange.m in data
        Sbar
    end
    methods
 
        function [a, F_ar, F_av] = compAccelPar(obj,X_state,~)
            time_et = obj.et;
            TEI = cspice_pxform( 'J2000', 'ITRF93', time_et);
            P_ecef = TEI*X_state(1:3);
            [g_ecef,Jac] = GravAccelPartialsExterior_mex(obj.deg, obj.ord, obj.Re/1000, obj.mu_e*1e-9, P_ecef(1:3)'/1000, obj.Cbar, obj.Sbar);
            a = TEI'*g_ecef*1e3;
            
            F_ar = TEI'*Jac*TEI;
            
            F_av = zeros(3,3); 
        end            
        
        
        %%%%%%%%%%%%%%% Acceleration function hasnt been modified to take
        %%%%%%%%%%%%%%% order as an input yet
        function a = compAccel(obj,X_state,~)            % X_state = [vel;pos]
            time_et = obj.et;
            TEI = cspice_pxform( 'J2000', 'ITRF93', time_et);
            P_ecef = TEI*X_state(1:3);
            g_ecef = GravAccelExterior_mex(obj.deg, obj.ord, obj.Re/1000, obj.mu_e*1e-9, P_ecef(1:3)'/1000, obj.Cbar, obj.Sbar);
            a = TEI'*g_ecef*1e3;
            
        end
%         
%         function F_ar = compParPos(obj,X_state,~)            % compute jacobian
% %             time_et = t+obj.et;
%             time_et = obj.et;
%             TEI = cspice_pxform( 'J2000', 'ITRF93', time_et);
%             P_ecef = TEI*X_state(1:3);
%             [~,Jac,~,~] = GravAccelPartialsExterior_mex(obj.deg, obj.Re/1000, obj.mu_e*1e-9, P_ecef(1:3)'/1000, obj.Cbar, obj.Sbar);
%             F_ar = TEI'*Jac*TEI;   
%         end
        

         
        function F_av = compParVel(~,~,~)
            F_av = zeros(3,3);
        end
    end
    
    
end