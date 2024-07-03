classdef SRP_Fourier_class
    properties
        Mass                              % Area to mass ratio of satellite
        AU                                % Astronomical Unit
        P_sun                             % Radiation pressure
        X_sun                             % Position of sun w.r.t Earth in J2000 frame
        Re                                % Radius of Earth
        Rsun                              % Radius of Sun                            
        Fx
        Fy
        Fz
        SBF2ECI
    end
    
    methods
          function obj =  SRP_Fourier_class(mass, AU, P_sun, X_sun, R, Rsun, Fx, Fy, Fz, SBF2ECI)
            obj.Mass      = mass;
            obj.AU        = AU;
            obj.P_sun     = P_sun;
            obj.X_sun     = X_sun;
            obj.Re        = R;
            obj.Rsun      = Rsun;
            obj.Fx        = Fx;
            obj.Fy        = Fy;
            obj.Fz        = Fz;
            obj.SBF2ECI   = SBF2ECI;
        end       
        
        function [a, F_ar, F_av, Fa_Cr] = compAccelPar(obj,X_state,~)       % compute partial from positions
%             global sp_angle
            X = X_state(1:3,1);               % X_state = [vel;pos] w.r.t to primary
            X_c = obj.X_sun;
            X_rel = X_c - X;                 % third body w.r.t spacecraft
            r_rel = norm(X_rel);
            Ps  = obj.P_sun;
            mass = obj.Mass;
            au = obj.AU;
            fs = shadow_func(obj, X_state);
            a = -obj.SBF2ECI*Ps*fs*au^2/(mass*r_rel^2)*[obj.Fx;obj.Fy;obj.Fz];
            
            F_ar = -Ps*fs*au^2*[obj.Fx;obj.Fy;obj.Fz]*(2*X_rel'/r_rel^4);
            F_av = zeros(3,3);
            Fa_Cr = -obj.SBF2ECI*Ps*fs*au^2/(mass*r_rel^2);
        end
        
        
        function fs = shadow_func(obj, X_state, ~)
            X = X_state(1:3,1);               % X_state = [vel;pos] w.r.t to primary
            X_c = obj.X_sun;
            X_rel = X - X_c;                  %  spacecraft w.r.t third body
            a = asin(obj.Rsun/norm(X_rel));
            b = asin(obj.Re/norm(X));
            c = acos(X'*(X_rel)/(norm(X)*norm(X_rel)));
            if c <= abs(a-b)
                fs = 0;
            elseif c >= a+b
                fs = 1;
            elseif c < a+b && c > abs(a-b)
                x = (c^2+a^2-b^2)/(2*c);
                y = sqrt(a^2-x^2);
                A = a^2*acos(x/a) + b^2*acos((c-x)/b) - c*y;
                fs = 1-A/(pi*a^2);
            else 
                fs = 0;
            end
        end
        
    end
    
    
end