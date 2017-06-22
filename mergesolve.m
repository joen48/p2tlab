function [X, output] = mergesolve(polysys, varargin)
%MERGESOLVE Recursive polynomial root-finding algorithm. 
%   [U, output] = mergesolve(polysys, varargin) recursively solves polysys, 
%   which is a cell containing the coefficients and monomial exponents 
%   of the set of polynomial equations. mergesolve(polysys) recursively divides
%   the null space of the Macaulay matrix in block terms and conquers the roots 
%   that are found from a CPD at the lowest level.
%
%   The structure output contains the output of the algorithm used for
%   computing the decomposition and additional information on the options used. 
%
%   Options are:
%      options.degree               If present, the degree of the Macaulay
%                                   matrix that is constructed to solve the
%                                   set of polynomial equations.
%      options.Nsol                 If present, the number of solutions,
%                                   e.g., when solving an over-constrained
%                                   set. This option will override
%                                   options.RecursiveOrthogonalization to false.
%      options.sparse               If true, compute the numerical basis for 
%                                   the null space of the Macaulay matrix 
%                                   using a sparse QR algorithm. Uses SuiteSparseQR.
%                                   False by default. This option will override
%                                   options.RecursiveOrthogonalization to false.
%      options.RecursiveOrthogonalization If true, update the numerical
%      = [true|false|'auto']        basis for the null space of the Macaulay
%                                   matrix using the recursive
%                                   orthogonalization theorem. If 'auto',
%                                   update recursively if the degree is
%                                   sufficiently high. By default, 'auto'
%                                   is used.
%      options.Compression          Compress the tensorized null space of
%      = ['svd'|'none']             the Macaulay matrix by means of an
%                                   (ML)SVD. By default, 'svd' is used.
%      options.Complex              If true, compute a complex solution.
%      = [true|false]
%      options.BTD                  The main algorithm for computing the BTD.
%      = [@btd_nls|@btd_minf]       By default, @btd_nls is used.
%      options.CPD                  The main algorithm for computing the CPD.
%      = [@cpd||@cpd_gevd|@cpd_nls] By default, @cpd is used.
%
%   See also:
%       poly_btd, poly_cpd, polysys2tens

global BTD CPD

p = inputParser;
p.addOptional('degree', []);
p.addOptional('Nsol', []);
p.addOptional('sparse', false);
p.addOptional('RecursiveOrthogonalization', false);
p.addOptional('Complex', true);
p.addOptional('Compression', 'svd');
p.addOptional('BTD', @btd_nls);
p.addOptional('CPD', @cpd);
p.addOptional('Tol', 1e-6);         % should actually depend on the accuracy of the polynomial coefficients
p.addOptional('MaxIter', 500);      % 200, 500, 10000
p.parse(varargin{:});

options = cell2struct(struct2cell(p.Results), fieldnames(p.Results));

% number of equations
s = size(polysys, 1);
% number of variables
n = size(polysys{1,2}, 2);
if (s < n)
    error('mergesolve:solution', 'The solution set is not 0-dimensional.')
end
% this vector will contain the degree of each equation
d_i = zeros(1,s);
for i = 1:s
    d_i(i) = max(sum(polysys{i,2}, 2)); 
end
if ~isempty(options.Nsol)
    m = options.Nsol;
    options.RecursiveOrthogonalization = false;
else
    % Bézout number
    m = prod(d_i);
    options.Nsol = m;   % update for output
end
dbound = sum(d_i)-n;
if isempty(options.degree)
    options.degree = dbound+1;
end

if options.sparse
    options.RecursiveOrthogonalization = false;
end

T = polysys2tens(polysys, 1, struct('degree', options.degree, 'm', m, ...
    'sparse', options.sparse, 'RecursiveOrthogonalization', options.RecursiveOrthogonalization, ...
    'd0', max(d_i), 'Compression', options.Compression));

% Recursively decompose T
if ~islogical(options.Complex)
    warning('polymerge:Complex', 'Invalid Complex option. Will use true instead.')
    options.Complex = true;
end
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if ~(xsfunc(options.BTD)&&xsfunc(options.CPD))
    error('polymerge:Algorithm', 'Invalid algorithm specified.')
else
    BTD = options.BTD;
    CPD = options.CPD;
end

[X, output] = div_con(T, randn(n+1,m), struct('Complex', options.Complex, ...
    'MaxIter', options.MaxIter, 'TolX', options.Tol, 'TolFun', options.Tol^2));
               
clear BTD CPD

output.Options = options;

end

function [X, output] = div_con(T, U1, options)
% Recurisively decompose T. The CPD is initialized with U1 in its first mode. 

global BTD CPD

% number of solutions
m = size(T,3);

if (m == 2)     % stop
    if strcmpi(func2str(CPD), 'cpd') || strcmpi(func2str(CPD), 'cpd_gevd')
        [U, output.CPD] = CPD(T, m, options, 'IsReal', ~options.Complex);
    else
        if options.Complex
            U0 = cpd_rnd(size(T), m, 'Real', @randn, 'Imag', @randn);
        else
            U0 = cpd_rnd(size(T), m, 'Real', @randn);
        end
        U0{1} = U1;
        [U, output.CPD] = CPD(T, U0, options, 'LargeScale', (sum(size(T))*m > 1e2));
    end
    U = rescale(U, options.TolX);
    X = U{1};
    return
elseif (m == 1)
    U = rescale({T(:,:,1), 1, 1}, options.TolX);
    X = U{1};
    return
end

% divide
l = floor(m/2);
r = m-l;
size_core = {repmat(l, [1 3]), repmat(r, [1 3])};
if options.Complex
    U0 = btd_rnd(size(T), size_core, 'Real', @randn, 'Imag', @randn);
else
    U0 = btd_rnd(size(T), size_core, 'Real', @randn);
end
U0{1}{1} = U1(:,1:l);
U0{2}{1} = U1(:,l+1:end);
[U, output.BTD] = BTD(T, U0, options);

% conquer
if (l == 1)
    U{1} = rescale(U{1}, options.TolX);
    Xl = U{1}{1};
else
    L = compress23(U{1});
    [Xl, output.left] = div_con(L, U{1}{1}, options);
end
if (r == 1)
    U{2} = rescale(U{2}, options.TolX);
    Xr = U{2}{1};
else
    R = compress23(U{2});
    [Xr, output.right] = div_con(R, U{2}{1}, options);
end
X = [Xl Xr];

end

function T = compress23(U)
% [U, ~, ~] = mlsvd(T);
[Q2, ~] = qr(U{2},0); 
U{2} = Q2'*U{2};
[Q3, ~] = qr(U{3},0);
U{3} = Q3'*U{3};
T = btdgen({U});
end


