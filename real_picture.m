function real_picture(polysys, rec)
%REALPICTURE Real picture of the zero set of a bivariate polynomial system.
%   realpicture(polysys, rec) draws the real picture of the zero set of
%   polysys, which is a 2 x 2 cell containing the coefficients and monomial 
%   exponents of a bivariate set of polynomial equations on the rectangular
%   grid rec = [xmin, xmax; ymin, ymax].
%
%   See also:
%       randpoly

% number of equations
s = size(polysys, 1);
if (s~=2)
    error('real_picture:bivariate', 'Only square bivariate polynomial systems are supported.'); 
end

% rectangular grid
x = linspace(rec(1,1), rec(1,2));
y = linspace(rec(2,1), rec(2,2));
[X, Y] = meshgrid(x, y);

% Color
c = {[0,0.447,0.741], [0.850,0.325,0.098]};

figure
for n = 1:s
    if ~(size(polysys{n,2},2)==2)
        error('real_picture:bivariate', 'Only square bivariate polynomial systems are supported.'); 
    end
    Z = zeros(size(X));
    for i = 1:size(Z,1)
        for j = 1:size(Z,2)
            Z(i,j) = sum(polysys{n,1} .* (X(i,j).^(polysys{n,2}(:,1))) .* (Y(i,j).^(polysys{n,2}(:,2))));
        end
    end
    contour(X, Y, Z, [1e-3 1e-3], 'Color', c{n}), hold on 
%     surf(X, Y, Z), hold on
end
hold off, grid on, box on, axis equal
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
