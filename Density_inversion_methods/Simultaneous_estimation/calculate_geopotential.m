function [V] = calculate_geopotential(deg_max, ord_max, a, mu, recef, C, S)

% INPUTS:
%   recef, size [3xN], Earth-centered Earth-fixed position vector for N data points
%   nmax,  size   [1], degree of truncation for gravitational spherical harmonic expansion
% OUTPUTS:
%   V,     size [1xN], Gravitational geopotential

%% Define constants
sqrt2 = sqrt(2);

%% calculate geocentric lat/lon/r
gc_r = sqrt(sum(recef.^2,1));
gc_lat = asind(recef(3,:)./gc_r);
gc_lon = atan2d(recef(2,:),recef(1,:));

%% Calculate gravitational potential
sinlat = sind(gc_lat);
[coslon,sinlon] = deal(zeros(deg_max,length(gc_lon)));
for m = 1:ord_max
	coslon(m,:) = cosd(m*gc_lon);
	sinlon(m,:) = sind(m*gc_lon);
end

% Accumulate (1+SUM(degree)(a/r)^n SUM(order) Cnm*Ynm(th,lam)) 
V = 1;
for n = 2:deg_max
	% The MATLAB legendre function provides:
	% legendre(n,x,'sch'); % schmidt normalization:
	% SP(N,M;X) = P(N,X), M = 0
        %           = (-1)^M * sqrt(2*(N-M)!/(N+M)!) * P(N,M;X), M > 0
	% or
	% legendre(n,x,'norm'); % 'full' normalization:
	% NP(N,M;X) = (-1)^M * sqrt((N+1/2)*(N-M)!/(N+M)!) * P(N,M;X)
	%
	% The desired normalization for geodesy:
	% GP(N,M;X) = sqrt(2*N+1)*P(N,X), M = 0
        %           = sqrt(2*(2*N+1)*(N-M)!/(N+M)!) * P(N,M;X), M > 0
	% so, the easiest conversion is the following:
	% GP(N,M;X) = sqrt(2) * NP(N,M;X), M = 0
	%           = 2 * NP(N,M;X) , M > 0
	% 
	GP = legendre(n,sinlat,'norm');

	% transform from tide-free to zero-tide (if desired??)
	%C(3,1) = C(3,1) - 3.11080e-8 * 0.3 / sqrt(5);

	% M = 0
	norm_factor = sqrt2;
	V = V + norm_factor*C(n+1,1)*(a./gc_r).^n .* GP(1,:);

	% M = 1:N
	for m = 1:n
		norm_factor = 2; % note: no Condon-Shortly phase here
		V = V + norm_factor*(a./gc_r).^n .* ( C(n+1,m+1)*coslon(m,:) + S(n+1,m+1)*sinlon(m,:) ) .* GP(m+1,:);
	end

end

V = mu*V./gc_r;

return
