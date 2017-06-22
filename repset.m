function over_polysys = repset(polysys, SNR, N, intrinsic)
%REPSET Over-constrained set of polynomial equations.
%   over_polysys = repset(polysys, SNR, N, instrinsic) creates an 
%   overdetermined set of polynomial equations consisting of N noisy copies 
%   of polysys. SNR is the Signal-to-Noise Ratio. If intrinsic is true, 
%   the coefficients of the intrinsic support are left untouched.
%
% See also:
%   poly_cpd

if (nargin < 3), N = 1; end
if (nargin < 4), intrinsic = true; end

% number of equations
s = size(polysys, 1);

over_polysys = cell(N*s,2);
for i = 1:s
    for j = 1:N
        if intrinsic
            over_polysys{(i-1)*N+j,1,1} = polysys{i,1} .* (abs(polysys{i,1})==1) ...
                + noisy(polysys{i,1}, SNR) .* (abs(polysys{i,1})~=1);
        else
            over_polysys{(i-1)*N+j,1} = noisy(polysys{i,1}, SNR);
        end
        over_polysys{(i-1)*N+j,2} = polysys{i,2};
    end
end

end


    