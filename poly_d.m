function polysys = poly_d(polysys, delta)
%POLY_D Partial derivative of a set of polynomial equations.
%   polysys = poly_d(polysys, delta) applies the differential functional 
%   delta to each equation in polysys, which is a cell containing the 
%   coefficients and monomial exponents of the set of polynomial equations. 

% number of equations
s = size(polysys, 1);
% number of variables
n = size(polysys{1,2}, 2);
if ~(size(delta, 2) == n)
    error('real_picture:delta', 'Size of the differential functional should match the number of (affine) variables.')
end

for i = 1:s
    p = size(polysys{i,1}, 1);
    for l = 1:p
        for j = 1:n
            for k = 1:delta(j)
                d = polysys{i,2}(l,j);
                polysys{i,1}(l) = d*polysys{i,1}(l);
                polysys{i,2}(l,j) = max(0, d-1);
            end
        end
    end
end
                

