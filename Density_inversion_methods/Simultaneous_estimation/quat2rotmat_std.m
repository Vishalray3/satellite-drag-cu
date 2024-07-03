function [rot_matrix] = quat2rotmat_std(Q)
%function quaternions takes a set of quaternions
%as input, and outputs a rotations matrix
% Convention below:
%q1 = scalar (convention is of the form 'VVLH -> SBF')
%q2 = x vector-component
%q3 = y vector-component
%q4 = z vector-component
% In order to convert, do this:


rot_matrix(1,1)= Q(2)^2-Q(3)^2-Q(4)^2+Q(1)^2;
rot_matrix(2,2)=-Q(2)^2+Q(3)^2-Q(4)^2+Q(1)^2;
rot_matrix(3,3)=-Q(2)^2-Q(3)^2+Q(4)^2+Q(1)^2;
rot_matrix(1,2)=2*(Q(2)*Q(3) + Q(1)*Q(4));
rot_matrix(1,3)=2*(Q(2)*Q(4) - Q(1)*Q(3));
rot_matrix(2,1)=2*(Q(3)*Q(2) - Q(1)*Q(4));
rot_matrix(2,3)=2*(Q(3)*Q(4) + Q(1)*Q(2));
rot_matrix(3,1)=2*(Q(4)*Q(2) + Q(1)*Q(3));
rot_matrix(3,2)=2*(Q(4)*Q(3) - Q(1)*Q(2));
