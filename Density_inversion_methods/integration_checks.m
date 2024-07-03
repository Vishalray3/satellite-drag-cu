%% Integral check
clc
clearvars
s = 8;
r = 0.08;
n = 5;

% fun = @(x,s,n)exp(-cos(x).^2.*s^2).*cos(x).*cos(n*x);
% int_num = integral(@(x)fun(x,s,n),0,2*pi);
% int_an = pi*exp(-s^2/2)*(besseli((n+1)/2, -s^2/2) + besseli((n-1)/2, -s^2/2) );
% int_num 
% 
% fun = @(x,s,n)(erf(cos(x)*s)).*cos(n*x);
% int_num = integral(@(x)fun(x,s,n),0,2*pi);
% int_an = 2*sqrt(pi)*s*exp(-s^2/2)*( (besseli((n-1)/2, -s^2/2) - besseli((n+1)/2, -s^2/2))/(n));
% int_num - int_an

% fun = @(x,s,n)(1+erf(cos(x)*s)).*cos(x).*cos(n*x);
% int_num = integral(@(x)fun(x,s,n),0,2*pi);
% int_an = sqrt(pi)*s*exp(-s^2/2)*( (besseli((n)/2, -s^2/2) - besseli((n+2)/2, -s^2/2))/(n+1) + (besseli((n-2)/2, -s^2/2) - besseli((n)/2, -s^2/2))/(n-1));
% int_num - int_an

fun = @(x,s,n)(1+erf(cos(x)*s)).*cos(x).^2.*cos(n*x);
int_num = integral(@(x)fun(x,s,n),0,2*pi);
int_an = sqrt(pi)*s*exp(-s^2/2)*( (besseli((n-1)/2, -s^2/2) - besseli((n+1)/2, -s^2/2))/(n) + (besseli((n+1)/2, -s^2/2) - besseli((n+3)/2, -s^2/2))/(2*(n+2))...
          + (besseli((n-3)/2, -s^2/2) - besseli((n-1)/2, -s^2/2))/(2*(n-2)));
int_num - int_an