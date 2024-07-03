classdef TwoBody
    properties
        mu_e                                % gravitational constant
        
    end
    
    methods

        function [a, F_ar, F_av] = compAccelPar(obj,X_state,~)   % compute partial from positions
            X = X_state(1:3);                 % X_state = [vel;pos]
            x = X(1);
            y = X(2);
            z = X(3);
            r = norm(X);
            mue = obj.mu_e;
            a(1,1) = -(mue*x/r^3);
            a(2,1) = -(mue*y/r^3);
            a(3,1) = -(mue*z/r^3);
            
            r_cap = X/r;
            F_ar = mue/(r)^3*( 3*(r_cap*r_cap') - eye(3));
            F_av = zeros(3,3);
        end
        
        function a = compAccel(obj,X_state,~)   % compute partial from positions
            X = X_state(1:3);                 % X_state = [vel;pos]
            x = X(1);
            y = X(2);
            z = X(3);
            r = norm(X);
            mue = obj.mu_e;
            a(1,1) = -(mue*x/r^3);
            a(2,1) = -(mue*y/r^3);
            a(3,1) = -(mue*z/r^3);
        end
        
        function F_ar = compParPos(obj,X_state,~)            % compute jacobian
            X = X_state(1:3);
            mue = obj.mu_e;
            r = norm(X);
            r_cap = X/r;
            F_ar = mue/(r)^3*( 3*(r_cap*r_cap') - eye(3));
            
        end
        function F_amu = compParMu(~,X_state,~)
            X = X_state(1:3);
            r = norm(X);
            F_amu = -X/r^3;
        end
        
        function F_av = compParVel(obj, X_state, ~)
            F_av = zeros(3,3)*X_state(1:3);
        end
    end
    
    
end