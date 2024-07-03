function [ b_bar_real , b_bar_imag ] = GetBnmNormalizedExterior(n_degree,R_ref,r_sat,x_sat,y_sat,z_sat)
% =============================================================================
% Author: Siamak Hesar             Date: 03/28/2015
% Advisor: Dr. Scheeres (The University of Colorado at Boulder)
%
% File name: GetBnmNormalizedExterior.m
%
% Description:
%       This function compute the basis functions for the exterior gravity
%       field. This is a test matlab code. The actual code is in C, mex.
%
% Inputs:
%	- n_degree,R_ref,r_sat,x_sat,y_sat,z_sat
%
% Outputs:
%	- b_bar_real,b_bar_imag
%
% References:
%	- 1: R. A. Werner, "Evaluating Descent and Ascent Trajectories Near Non-Spherical Bodies", Technical Support Package		
%
% Mathematical Formulation:
%	-
% Note:
%  - 
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
% - 03/28/2015		Siamak Hesar		original version
%	
%% ============================================================================  
% -- Pre-allocation -- 
b_bar_real = zeros(n_degree);
b_bar_imag = zeros(n_degree);


%% Vertical Recurrences -- 

for mm = 1:n_degree
    
    m = mm - 1;

    for nn = mm:n_degree
        
        n = nn - 1;

        if mm == nn

            if m == 0

                b_bar_real(1,1) = R_ref/r_sat;
                b_bar_imag(1,1) = 0.0;

            else

                if (n == 1)

                    delta_1_n = 1.0;

                else

                    delta_1_n = 0.0;

                end

                b_bar_real(nn,nn) = sqrt( (1.0 + delta_1_n)*(2.0*n + 1.0)/(2.0*n) ) * (R_ref/r_sat) * ( x_sat/r_sat*b_bar_real(nn-1,nn-1) - y_sat/r_sat*b_bar_imag(nn-1,nn-1) );
                b_bar_imag(nn,nn) = sqrt( (1.0 + delta_1_n)*(2.0*n + 1.0)/(2.0*n) ) * (R_ref/r_sat) * ( y_sat/r_sat*b_bar_real(nn-1,nn-1) + x_sat/r_sat*b_bar_imag(nn-1,nn-1) );

            end

        else

            if ( n >= 2 )

                b_bar_real(nn,mm) = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(R_ref/r_sat)*(z_sat/r_sat)*b_bar_real(nn-1,mm) - sqrt( (2.0*n + 1.0)*( (n - 1.0)*(n - 1.0) - m*m )/( (2.0*n - 3.0)*(n*n - m*m) ) )*(R_ref/r_sat)*(R_ref/r_sat)*b_bar_real(nn-2,mm);
                b_bar_imag(nn,mm) = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(R_ref/r_sat)*(z_sat/r_sat)*b_bar_imag(nn-1,mm) - sqrt( (2.0*n + 1.0)*( (n - 1.0)*(n - 1.0) - m*m )/( (2.0*n - 3.0)*(n*n - m*m) ) )*(R_ref/r_sat)*(R_ref/r_sat)*b_bar_imag(nn-2,mm);

            else

                b_bar_real(nn,mm) = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(R_ref/r_sat)*(z_sat/r_sat)*b_bar_real(nn-1,mm);
                b_bar_imag(nn,mm) = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(R_ref/r_sat)*(z_sat/r_sat)*b_bar_imag(nn-1,mm);

            end

        end

    end

end
