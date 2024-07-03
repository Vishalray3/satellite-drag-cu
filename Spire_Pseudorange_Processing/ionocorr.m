function [PRIF, iono] = ionocorr (C1, f1, P2, f2)
%==================================================
%   Author: Suood Alnaqbi 
%   10/11/2020 
%==================================================
% This function calculates the ionospheric correction using dual-frequency
% approach. 
%==================================================
%   Inputs: 
%       C1           psuedorange measruments using L1 frequency 
%       f1           L1 carrier frequency in Hz 
%       P2           Psuedorange measurments using L2 frequency 
%       f2           L2 carrier frequency in Hz 
%   Outputs
%       PRIF         Ionosphere-free range 
%       iono
%==================================================


PRIF=(f1^2*C1-f2^2*P2)/(f1^2-f2^2); % Ionosphere-free range in m 
TEC=(f1^2*f2^2)/(40.3*(f1^2-f2^2))*(P2-C1); % Total electron content, 10^16 electrons /m^2 
Ip1=40.3*TEC/f1^2; % Ionospheric correction for L1 measurments
Ip2=40.3*TEC/f2^2; % Ionospheric correction for L2 measurments 
iono=[Ip1;Ip2]; % Ionospheric corrections for both L1 and L2 












