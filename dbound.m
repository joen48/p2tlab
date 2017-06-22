function d = dbound(polysys)
%DBOUND Get upper bound for the degree of regularity.
%   d = dbound(polysys) returns the upper bound for the degree of
%   regularity of polysys, which is a cell containing the coefficients and
%   monomial exponents of a set of polynomial equations. 
%
% See also:
%   poly_cpd, polysys2tens

% number of equations
s = size(polysys, 1);
% number of variables
n = size(polysys{1,2}, 2);
d = 0;
for i = 1:s
    d = d + max(sum(polysys{i,2}, 2));    % d* = d* + d_i 
end
d = d-n;
