function m = bezout(polysys)
%BEZOUT Bézout number.
%   m = bezout(polysys) returns the Bézout number of polysys, 
%   which is a cell containing the coefficients and monomial exponents of
%   a set of polynomial equations.
%
% See also:
%   poly_cpd, polysys2tens

% number of equations
s = size(polysys, 1);
% number of variables
n = size(polysys{1,2}, 2);
m = 1;
for i = 1:s
    m = m * max(sum(polysys{i,2}, 2));    % m = m * d_i 
end
