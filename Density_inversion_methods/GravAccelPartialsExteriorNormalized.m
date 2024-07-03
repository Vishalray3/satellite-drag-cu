function [ Accel, Jacobian, Partial_C, Partial_S] = GravAccelPartialsExteriorNormalized(n_degree, R_ref, mu, r_vec, Cbar, Sbar)
%%=========================================================================
% File name: GravAccelPartialsExteriorNormalized.m
%
% Author: Siamak Hesar             3/22/2015
%
% This function is adapted from the original code written by Yu Takahashi 
% called AccelInteriorPotential_mex.c. However I implemented Hotine's approach
% in computing first and second partial derivatives of the potential w.r.t
% the Cartesian coordinates. The original code used Cunningham's method.
% 
% Edit by Vishal Ray               06/09/2019   
% : Added the input of gravity field order alongwith the degree and
% modified the code accordingly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%
% Description:
%
%  This function computes the acceleration, and dynamics 
%  matrix for the exterior potential. It returns the full partial matrix 
%  of acceleration w.r.t state and gravity harmonics partials.
%
% Inputs:
%
%     n_degree   (n.d.)     : Degree of the spherical harmonics
%
%     R_ref      (km)       : Reference distance. usually the reference radius
%
%     mu         (km^3/sec^2): Gravitational Parameter of the asteroid
%
%     r_vec      (km)       : Field point (spacecraft) vector = (x_sat, y_sat, z_sat)
%
%     Cbar       (n.d.)     : Normalized C spherical harmonics
%
%     Sbar       (n.d.)     : Normalized S spherical harmonics
%
% Outputs:
%
%     Accel              : Acceleration by basis \bar{K}_{nm}^i (normalized)
%
%     Jacobian           : STM by \bar{b}_{nm}^i
%
%     Partial_C          : Partial Accel Partial C
%
%     Partial_S          : Partial Accel Partial S
%
% Assumptions/References:
%	- 1: S. V. Bettadpur, "Hotine's geopotential formulation: revisited", Bulletin Geodesique (1995) 69:i35-142
%   - 2: R. A. Werner, "Evaluating Descent and Ascent Trajectories Near Non-Spherical Bodies", Technical Support Package
%	- 3: L. E. Cunningham, "On the computation of the spherical harmonic terms needed during the numerical integration of the orbital motion of an artificial satellite"
%	
%  - Mathematical Formulation
%
%    ~ Note the following definitions
%
%       (1) b_{n,m}^e       = (R_ref/r)^{n+1} * Pnm * (cos(m*lambda); sin(m*lambda))
%
%       (2) \bar{b}_{n,m}^e = \sqrt{(2 - \delta_{0,m}) * (2n + 1) * (n - m)! / (n + m)!} * (R_ref/r)^{n+1} * Pnm * (cos(m*lambda); sin(m*lambda))
%
%       (3) c_{n,m}^e       = ( 2 - \delta_{0,m} )*(n - m)!/(n + m)!*(r'/R_ref)^n * Pnm * (cos(m*lambda'); sin(m*lambda'))
%
%       (4) \bar{c}_{n,m}^e = \sqrt{ (2 - \delta_{0,m})*(n - m)!/( (2n + 1)*(n + m)! ) } *(r'/R_ref)^n * Pnm * (cos(m*lambda'); sin(m*lambda'))
%
%      where ' indicates the parameters of the differential mass.
%
%    ~ Note that these expressions can be considered as imaginary numbers. That is,
%
%       (1) b_{n,m}^i       = (R_ref/r)^{n+1} * Pnm * e^{i*m*lambda}
%
%       (2) \bar{b}_{n,m}^i = \sqrt{(2 - \delta_{0,m}) * (2n + 1) * (n - m)! / (n + m)!} * (R_ref/r)^{n+1} * Pnm * e^{i*m*lambda}
%
%       (3) c_{n,m}^i       = ( 2 - \delta_{0,m} )*(n - m)!/(n + m)! * (r'/R_ref)^n * Pnm * e^{i*m*lambda'}
%
%       (4) \bar{c}_{n,m}^i = \sqrt{ (2 - \delta_{0,m})*(n - m)!/( (2n + 1)*(n + m)! ) } * (r'/R_ref)^n * Pnm * e^{i*m*lambda'}
%
%    ~ The normalization factor is defined as N, where
%
%        \bar{P}_{n,m} = N*P_{n,m}
%
%      Thus, N is immediately recognized as
%
%         N        = \sqrt{(2 - \delta_{0,m}) * (2n + 1) * (n - m)! / (n + m)!}
%
%    ~ The external potetial of the field point is computed as
%
%        (1) U^e       = \frac{G*M_ref*}{R_ref} * \sum^{\infty}_{n = 0} \sum^n_{m = 0} b_{n,m}^i * (1/M_ref) * \int_M c_{n,m}^i dm'
%
%        (2) U^e       = \frac{G*M_ref*}{R_ref} * \sum^{\infty}_{n = 0} \sum^n_{m = 0} \bar{b}_{n,m}^i * (1/M_ref) * \int_M \bar{c}_{n,m}^i dm'
%
%    ~ This function has the following recursive formulae:
%
%      -------------------------------------------------------------------
%
%      (1) Basis function          : b_{0,0}^e   = (R_ref/r) * (1; 0)
%
%      (2) Diagonal recurrences    : b_{n,n}^e   = (2*n - 1) * (R_ref/r) * (x/r, -y/r; y/r, x/r)*b_{n-1,n-1}^e
%
%      (3) Subdiagonal recurrences : b_{n,n-1}^e = (2*n - 1) * (R_ref/r) * (z/r) * b_{n-1,n-1}^e
%
%      (4) Vertical recurrences    : b_{n,m}^e   = \frac{2*n - 1}{n - m}* (R_ref/r) * (z/r) * b_{n-1,m}^e - \frac{n + m - 1}{n - m}*(R_ref/r)^2*b_{n-2,m}^e
%      -------------------------------------------------------------------
%
%      (1) Basis function          : \bar{b}_{0,0}^e   = (R_ref/r) * (1; 0)
%
%      (2) Diagonal recurrences    : \bar{b}_{n,n}^e   = \sqrt{(1 + \delta_{1,n})*(2n + 1)/(2n)}* (R_ref/r) * (x/r, -y/r; y/r, x/r)*\bar{b}_{n-1,n-1}^e
%
%      (3) Subdiagonal recurrences : \bar{b}_{n,n-1}^e = \sqrt{2*n - 1}* (R_ref/r) * (z/R_ref)*\bar{b}_{n-1,n-1}^e
%
%      (4) Vertical recurrences    : \bar{b}_{n,m}^e   = \frac{4n^2 - 1}{n^2 - m^2}* (R_ref/r) * (z/R_ref)*\bar{b}_{n-1,m}^e - \frac{(2n + 1)*( (n - 1)^2 - m^2 )}{(2n - 3)*(n^2 - m^2)}*(R_ref/r)^2*\bar{b}_{n-2,m}^e
%
%      -------------------------------------------------------------------
%
%      (1) Basis function          : c_{0,0}^e   = (1; 0)
%
%      (2) Diagonal recurrences    : c_{n,n}^e   = (1 + \delta_{1,n})/(2n)*(x'/R_ref, -y'/R_ref; y'/R_ref, x'/R_ref)*c_{n-1,n-1}^e
%
%      (3) Subdiagonal recurrences : c_{n,n-1}^e = (z'/R_ref)*c_{n-1,n-1}^e
%
%      (4) Vertical recurrences    : c_{n,m}^e   = \frac{2n - 1}{n + m}*(z'/R_ref)*c_{n-1,m}^e - \frac{n - m - 1}{n + m}*(r'/R_ref)^2*c_{n-2,m}^e
%
%      -------------------------------------------------------------------
%
%      (1) Basis function          : \bar{c}_{0,0}^e   = (1; 0)
%
%      (2) Diagonal recurrences    : \bar{c}_{n,n}^e   = (2n - 1)*\sqrt{ (1 + \delta_{1,n})/( (2n)*(2n + 1) ) }*(x'/R_ref, -y'/R_ref; y'/R_ref, x'/R_ref)*\bar{c}_{n-1,n-1}^e
%
%      (3) Subdiagonal recurrences : \bar{c}_{n,n-1}^e = \frac{(2n - 1)}{\sqrt{2n + 1}}*(z'/R_ref)*\bar{c}_{n-1,n-1}^e
%
%      (4) Vertical recurrences    : \bar{c}_{n,m}^e   = (2n - 1)*\sqrt{\frac{2n - 1}{(2n + 1)*(n^2 - m^2)}}*(z'/R_ref)*\bar{c}_{n-1,m}^e - \sqrt{\frac{(2n - 3)*( (n - 1)^2 - m^2 )}{(2n + 1)*(n^2 - m^2)}}*(r'/R_ref)^2*\bar{c}_{n-2,m}^e
%
%      -------------------------------------------------------------------
%
%      where subdiagonal recurrences are determined by m = n - 1 through the vertical recurrences.
%
%      -------------------------------------------------------------------
%
% Note:
%
%  - Modified from the original version 
%
% Dependencies:
%
%  - None
%
% Call
%
%  - None
%
% Called by
%
% - TBD
%
% Modification History:
%
%  27Feb11   Yu Takahashi   original version of AccelInteriorPotential_mex.c
%  
%  3/22/2015   Siamak Hesar    Mdified from the original code written by Yu Takahashi (AccelInteriorPotential_mex.c)
%                              to compute the normalized accelerations and full partials matrix for an exterior gravity field.
% 
%  3/22/2015   Siamak Hesar    Added statements for validating the mex function input types and sizes.
%% ========================================================================
% - Satellite position
    
x_sat = r_vec(1); y_sat = r_vec(2); z_sat = r_vec(3); % (km) x, y, z position vector of the spacecraft
r_sat = sqrt(x_sat*x_sat + y_sat*y_sat + z_sat*z_sat);            % (km) Norm of the position vector

% - Basis functions
[b_bar_real,b_bar_imag] = GetBnmNormalizedExterior(n_degree+3 , R_ref, r_sat, x_sat, y_sat, z_sat);

% - Pre-allocation
x_ddot = 0.0;       y_ddot = 0.0;       z_ddot = 0.0; 
ddU_dxdx = 0.0;     ddU_dxdy = 0.0;     ddU_dxdz = 0.0;
ddU_dydy = 0.0;     ddU_dydz = 0.0;     ddU_dzdz = 0.0;

%% First Partials of the Potential
 
K0 = 0.5*mu/R_ref^2;

for nn = 1 : n_degree + 1
    
    n = nn - 1;
%        if nn <= m_order + 1
%            kk = nn;
%        else
%            kk = m_order + 1;
%        end
    for mm = 1 : nn
        
        m = mm - 1;
        
        if m == 1
            delta_1_m = 1;
        else
            delta_1_m = 0;
        end
        
        K1 = sqrt((n+2)*(n+1)*(2*n+1)/2/(2*n+3));
        K2 = sqrt((n+m+2)*(n+m+1)*(2*n+1)/(2*n+3));
        K3 = sqrt(2*(n-m+2)*(n-m+1)*(2*n+1)/(2 - delta_1_m)/(2*n+3));
        
        if m == 0
            
            x_ddot = x_ddot - 2*K0 * ( Cbar(nn,mm)*K1*b_bar_real(nn+1,mm+1) );
            y_ddot = y_ddot - 2*K0 * ( Cbar(nn,mm)*K1*b_bar_imag(nn+1,mm+1) );
            z_ddot = z_ddot - 2*K0 * ( Cbar(nn,mm)*sqrt((n-m+1)*(n+m+1)*(2*n+1)/(2*n+3))*b_bar_real(nn+1,mm) );
            
        else
            
            x_ddot = x_ddot + K0 * ( -Cbar(nn,mm)*K2*b_bar_real(nn+1,mm+1) -Sbar(nn,mm)*K2*b_bar_imag(nn+1,mm+1) +Cbar(nn,mm)*K3*b_bar_real(nn+1,mm-1) +Sbar(nn,mm)*K3*b_bar_imag(nn+1,mm-1));
            y_ddot = y_ddot + K0 * ( -Cbar(nn,mm)*K2*b_bar_imag(nn+1,mm+1) +Sbar(nn,mm)*K2*b_bar_real(nn+1,mm+1) -Cbar(nn,mm)*K3*b_bar_imag(nn+1,mm-1) +Sbar(nn,mm)*K3*b_bar_real(nn+1,mm-1));
            z_ddot = z_ddot - 2*K0 * ( Cbar(nn,mm)*sqrt((n-m+1)*(n+m+1)*(2*n+1)/(2*n+3))*b_bar_real(nn+1,mm) +Sbar(nn,mm)*sqrt((n-m+1)*(n+m+1)*(2*n+1)/(2*n+3))*b_bar_imag(nn+1,mm) );
            
        end
        
    end
    
end

%% Second Partials of the Potential

K0 = 0.25*mu/R_ref^3;

for nn = 1 : n_degree + 1
    
    n = nn - 1;
    
%        if nn <= m_order + 1
%            kk = nn;
%        else
%            kk = m_order + 1;
%        end
    for mm = 1 : nn
        
        m = mm - 1;
        
        if m == 1
            delta_1_m = 1;
        else
            delta_1_m = 0;
        end
        
        if m == 2
            delta_2_m = 1;
        else
            delta_2_m = 0;
        end        
        
        % - Common multipliers
        K1 = sqrt((n+m+4)*(n+m+3)*(n+m+2)*(n+m+1)*(2*n+1)/(2*n+5));
        K2 = sqrt((n-m+2)*(n-m+1)*(n+m+2)*(n+m+1)*(2*n+1)/(2*n+5));
        K3 = sqrt(2*(n-m+4)*(n-m+3)*(n-m+2)*(n-m+1)*(2*n+1)/(2 - delta_2_m)/(2*n+5));
        K4 = sqrt((n+5)*(n+4)*(n+3)*(n+2)*(2*n+1)/(2*n+5));
        K5 = sqrt((n+3)*(n+2)*(n+1)* n *(2*n+1)/(2*n+5));
        K6 = sqrt((n+4)*(n+3)*(n+2)*(n+1)*(2*n+1)/2/(2*n+5));
        K7 = sqrt((2*n+1)/(2*n+5));
        K8 = sqrt((n-m+1)*(n+m+3)*(n+m+2)*(n+m+1)*(2*n+1)/(2*n+5));
        K9 = sqrt(2*(n+m+1)*(n-m+3)*(n-m+2)*(n-m+1)*(2*n+1)/(2 - delta_1_m)/(2*n+5));
        K10= sqrt((n+3)*(n+2)*(n+1)*(n+1)*(2*n+1)/2/(2*n+5));
        
        % Partial expressions
        if m == 0
            
            ddU_dxdx = ddU_dxdx + 2*K0 * ( Cbar(nn,mm)*K6*b_bar_real(nn+2,mm+2) -(n+2)*(n+1)*K7*Cbar(nn,mm)*b_bar_real(nn+2,mm) );
            ddU_dydy = ddU_dydy - 2*K0 * ( Cbar(nn,mm)*K6*b_bar_real(nn+2,mm+2) +(n+2)*(n+1)*K7*Cbar(nn,mm)*b_bar_real(nn+2,mm) );
            ddU_dxdy = ddU_dxdy + 2*K0 * ( Cbar(nn,mm)*K6*b_bar_imag(nn+2,mm+2) );
            ddU_dxdz = ddU_dxdz + 4*K0 * ( Cbar(nn,mm)*K10*b_bar_real(nn+2,mm+1) );
            ddU_dydz = ddU_dydz + 4*K0 * ( Cbar(nn,mm)*K10*b_bar_imag(nn+2,mm+1) );
            ddU_dzdz = ddU_dzdz + 4*K0 * ( Cbar(nn,mm)*K2 *b_bar_real(nn+2,mm) );
            
        elseif m == 1
            
            ddU_dxdx = ddU_dxdx +   K0 * ( Cbar(nn,mm)*K4*b_bar_real(nn+2,mm+2) +Sbar(nn,mm)*K4*b_bar_imag(nn+2,mm+2) -3*K5*Cbar(nn,mm)*b_bar_real(nn+2,mm) -  K5*Sbar(nn,mm)*b_bar_imag(nn+2,mm));
            ddU_dydy = ddU_dydy -   K0 * ( Cbar(nn,mm)*K4*b_bar_real(nn+2,mm+2) +Sbar(nn,mm)*K4*b_bar_imag(nn+2,mm+2) +  K5*Cbar(nn,mm)*b_bar_real(nn+2,mm) +3*K5*Sbar(nn,mm)*b_bar_imag(nn+2,mm));
            ddU_dxdy = ddU_dxdy -   K0 * ( Sbar(nn,mm)*K4*b_bar_real(nn+2,mm+2) -Cbar(nn,mm)*K4*b_bar_imag(nn+2,mm+2) +  K5*Sbar(nn,mm)*b_bar_real(nn+2,mm) +  K5*Cbar(nn,mm)*b_bar_imag(nn+2,mm));
            ddU_dxdz = ddU_dxdz + 2*K0 * ( Cbar(nn,mm)*K8*b_bar_real(nn+2,mm+1) +Sbar(nn,mm)*K8*b_bar_imag(nn+2,mm+1) -  K9*Cbar(nn,mm)*b_bar_real(nn+2,mm-1) -K9*Sbar(nn,mm)*b_bar_imag(nn+2,mm-1) );
            ddU_dydz = ddU_dydz - 2*K0 * ( Sbar(nn,mm)*K8*b_bar_real(nn+2,mm+1) -Cbar(nn,mm)*K8*b_bar_imag(nn+2,mm+1) +  K9*Sbar(nn,mm)*b_bar_real(nn+2,mm-1) -K9*Cbar(nn,mm)*b_bar_imag(nn+2,mm-1) );
            ddU_dzdz = ddU_dzdz + 4*K0 * ( Cbar(nn,mm)*K2*b_bar_real(nn+2,mm)   +Sbar(nn,mm)*K2*b_bar_imag(nn+2,mm));
            
        else
            
            ddU_dxdx = ddU_dxdx + K0 * (  Cbar(nn,mm)*K1*b_bar_real(nn+2,mm+2) +Sbar(nn,mm)*K1*b_bar_imag(nn+2,mm+2) -2*K2*Cbar(nn,mm)*b_bar_real(nn+2,mm) -2*K2*Sbar(nn,mm)*b_bar_imag(nn+2,mm) +K3*Cbar(nn,mm)*b_bar_real(nn+2,mm-2) +K3*Sbar(nn,mm)*b_bar_imag(nn+2,mm-2));
            ddU_dydy = ddU_dydy - K0 * (  Cbar(nn,mm)*K1*b_bar_real(nn+2,mm+2) +Sbar(nn,mm)*K1*b_bar_imag(nn+2,mm+2) +2*K2*Cbar(nn,mm)*b_bar_real(nn+2,mm) +2*K2*Sbar(nn,mm)*b_bar_imag(nn+2,mm) +K3*Cbar(nn,mm)*b_bar_real(nn+2,mm-2) +K3*Sbar(nn,mm)*b_bar_imag(nn+2,mm-2));
            ddU_dxdy = ddU_dxdy + K0 * ( -Sbar(nn,mm)*K1*b_bar_real(nn+2,mm+2) +Cbar(nn,mm)*K1*b_bar_imag(nn+2,mm+2) +K3*Sbar(nn,mm)*b_bar_real(nn+2,mm-2) -K3*Cbar(nn,mm)*b_bar_imag(nn+2,mm-2));
            ddU_dxdz = ddU_dxdz + 2*K0 * ( Cbar(nn,mm)*K8*b_bar_real(nn+2,mm+1) +Sbar(nn,mm)*K8*b_bar_real(nn+2,mm+1) -  K9*Cbar(nn,mm)*b_bar_real(nn+2,mm-1) -K9*Sbar(nn,mm)*b_bar_imag(nn+2,mm-1) );
            ddU_dydz = ddU_dydz - 2*K0 * ( Sbar(nn,mm)*K8*b_bar_real(nn+2,mm+1) -Cbar(nn,mm)*K8*b_bar_imag(nn+2,mm+1) +  K9*Sbar(nn,mm)*b_bar_real(nn+2,mm-1) -K9*Cbar(nn,mm)*b_bar_imag(nn+2,mm-1) );
            ddU_dzdz = ddU_dzdz + 4*K0 * ( Cbar(nn,mm)*K2*b_bar_real(nn+2,mm)   +Sbar(nn,mm)*K2*b_bar_imag(nn+2,mm));
            
        end
               
    end
    
end

%% Partial Accel partial C&S

Ccounter = 1;
Scounter = 1;

K0 = 0.5*mu/R_ref^2;

Cno = (n+2)*(n+1) / 2;
Sno = n*(n+1) / 2;
Partial_C = zeros(3,Cno); 
Partial_S = zeros(3,Sno);

for nn = 1 : n_degree + 1
    
    n = nn - 1;
%        if nn <= m_order + 1
%            kk = nn;
%        else
%            kk = m_order + 1;
%        end
    for mm = 1 : nn
        
        m = mm - 1;
        
        if m == 1
            delta_1_m = 1;
        else
            delta_1_m = 0;
        end
        
        K1 = sqrt((n+2)*(n+1)*(2*n+1)/2/(2*n+3));
        K2 = sqrt((n+m+2)*(n+m+1)*(2*n+1)/(2*n+3));
        K3 = sqrt(2*(n-m+2)*(n-m+1)*(2*n+1)/(2 - delta_1_m)/(2*n+3));
        
        if m == 0
            
            Partial_C(1,Ccounter) = - 2*K0 * ( K1*b_bar_real(nn+1,mm+1) );
            Partial_C(2,Ccounter) = - 2*K0 * ( K1*b_bar_imag(nn+1,mm+1) );
            Partial_C(3,Ccounter) = - 2*K0 * ( sqrt((n-m+1)*(n+m+1)*(2*n+1)/(2*n+3))*b_bar_real(nn+1,mm) );
            Ccounter = Ccounter + 1;
                        
        else
            
            Partial_C(1,Ccounter) = K0 * ( -K2*b_bar_real(nn+1,mm+1) +K3*b_bar_real(nn+1,mm-1) );
            Partial_C(2,Ccounter) = K0 * ( -K2*b_bar_imag(nn+1,mm+1) -K3*b_bar_imag(nn+1,mm-1) );
            Partial_C(3,Ccounter) = - 2*K0 * ( sqrt((n-m+1)*(n+m+1)*(2*n+1)/(2*n+3))*b_bar_real(nn+1,mm) );
            Ccounter = Ccounter + 1;
            
            Partial_S(1,Scounter) = K0 * ( -K2*b_bar_imag(nn+1,mm+1) +K3*b_bar_imag(nn+1,mm-1));
            Partial_S(2,Scounter) = K0 * (  K2*b_bar_real(nn+1,mm+1) +K3*b_bar_real(nn+1,mm-1));
            Partial_S(3,Scounter) = -2*K0 * ( sqrt((n-m+1)*(n+m+1)*(2*n+1)/(2*n+3))*b_bar_imag(nn+1,mm) );
            Scounter = Scounter + 1;
                       
        end
                
    end
    
end

%% Arranging Outputs
Accel = [x_ddot;y_ddot;z_ddot];

Jacobian   = [ddU_dxdx, ddU_dxdy, ddU_dxdz;
              ddU_dxdy, ddU_dydy, ddU_dydz;
              ddU_dxdz, ddU_dydz, ddU_dzdz];

%% END OF FUNCTION
