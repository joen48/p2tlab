function [U, output] = poly_btd(polysys, size_core, varargin)
%POLY_BTD Multivariate polynomial subspace-finding by a BTD. 
%   [U, output] = poly_btd(polysys, size_core, varargin) decomposes the null
%   space of the Macaulay matrix of polysys, which is a cell containing the 
%   coefficients and monomial exponents of the set of polynomial equations. 
%   Each term in U{k} is a matrix of size (n+1) x size_core(k) that represents
%   the root x^{(k)} with multiplicity size_core(k).
%
%   poly_btd(polysys, G, varargin) initializes the core tensors in the BTD
%   U0{k}{4} with G{k}.
%
%   The structure output contains the output of the algorithm used for
%   computing the decomposition and additional information on the options used. 
%
%   Options are:
%      options.degree               If present, the degree of the Macaulay
%                                   matrix that is constructed to solve the
%                                   set of polynomial equations.
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
%      options.Algorithm            The main algorithm for computing the BTD.
%      = [@btd_nls|@btd_minf]
%
%   See also:
%       poly_cpd

p = inputParser;
p.addOptional('degree', []);
p.addOptional('sparse', false);
p.addOptional('RecursiveOrthogonalization', false);
p.addOptional('Complex', true);
p.addOptional('Compression', 'svd');
p.addOptional('Algorithm', @btd_nls);
p.addOptional('Tol', 1e-6);         % should actually depend on the accuracy of the polynomial coefficients
p.addOptional('MaxIter', 500);      % 200, 500, 10000
p.parse(varargin{:});

options = cell2struct(struct2cell(p.Results), fieldnames(p.Results));

% number of equations
s = size(polysys, 1);
% number of variables
n = size(polysys{1,2}, 2);
if (s < n)
    error('poly_btd:solution', 'The solution set is not 0-dimensional.')
end
% this vector will contain the degree of each equation
d_i = zeros(1,s);
for i = 1:s
    d_i(i) = max(sum(polysys{i,2}, 2)); 
end
% Bézout number
m = prod(d_i);
% upper bound for the degree of regularity
% dbound = sum(d_i)-n;

if ischar(options.RecursiveOrthogonalization) && strcmpi(options.RecursiveOrthogonalization, 'auto')
    % an estimated gain of d^3/n^3 operations when using the recursive orthogonalization scheme
    if (dbound^3/n^3 > 1)
        options.RecursiveOrthogonalization = true;  
    else
        options.RecursiveOrthogonalization = false;
    end
end
if options.sparse
    options.RecursiveOrthogonalization = false;
end
if ~islogical(options.RecursiveOrthogonalization)
    warning('poly_sd:RecursiveOrthogonalization', 'Invalid RecursiveOrthogonalization option. Will use false instead.')
    options.RecursiveOrthogonalization = false;
end
if ~ischar(options.Compression)
    warning('poly_cpd:Compression', 'Invalid Compression option. Will use ''svd'' instead.')
    options.Compression = 'svd';
end

if iscell(size_core)
    G = size_core;
    size_core = cellfun(@size, G, 'UniformOutput', false);
    for r = 1:numel(size_core)
        size_core{r} = [size_core{r} ones(1,3-numel(size_core{r}))];
    end
    L = max(cell2mat(size_core))-1;
elseif isnumeric(size_core)
    % From multiplicities to new size_core
    L = max(size_core)-1;
    mu = size_core;
    size_core = cell(size(mu));
    for r = 1:numel(size_core)
        size_core{r} = repmat(mu(r), [1 3]);
    end
else
    error('poly_btd:size_core', 'Invalid size_core.')
end

T = polysys2tens(polysys, L, struct('degree', options.degree, 'm', m, ...
    'sparse', options.sparse, 'RecursiveOrthogonalization', options.RecursiveOrthogonalization, ...
    'd0', max(d_i), 'Compression', options.Compression)); % size(Y,2) == m

if ~islogical(options.Complex)
    warning('poly_btd:Complex', 'Invalid Complex option. Will use true instead.')
    options.Complex = true;
end
if options.Complex
    U0 = btd_rnd(size(T), size_core, 'Real', @randn, 'Imag', @randn);
else
    U0 = btd_rnd(size(T), size_core, 'Real', @randn);
end   

if exist('G', 'var')
    % SDF
    model = struct;
    model.factorizations.poly_btd.data = T;
    model.factorizations.poly_btd.btd = cell([1 numel(G)]);
    for r = 1:numel(G)
        for N = 1:3
            model.variables.(sprintf('u%d%d', r, N)) = U0{r}{N};
            model.factors.(sprintf('U%d%d', r, N)) = sprintf('u%d%d', r, N);
        end
%         model.variables.(sprintf('g%d', r)) = U0{r}{4}(:,:,2:end);
%         model.factors.(sprintf('G%d', r)) = reshape({eye(size(U0{r}{4}, 1)); sprintf('g%d', r)}, 1, 1, 2);
        % Avoid references to variables
        model.factors.(sprintf('G%d', r)) = {G{r}, '(constant)'};
        model.factorizations.poly_btd.btd{r} = { sprintf('U%d%d', r, 1), ...
                                                 sprintf('U%d%d', r, 2), ...
                                                 sprintf('U%d%d', r, 3), ...
                                                 sprintf('G%d', r) };
    end
    % Investigate the resulting model
    sdf_check(model, 'print')
    % Solve it
    [sol, output.AlgorithmOutput] = sdf_nls(model, 'MaxIter', options.MaxIter, ...
        'TolX', options.Tol, 'TolFun', options.Tol^2, 'Complex', options.Complex);
    U = sol.factors;
else
    % Decomposition
    isfunc = @(f)isa(f,'function_handle');
    xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
    if ~xsfunc(options.Algorithm)
        error('poly_btd:Algorithm', 'Invalid Algorithm option.')
    end
    [U, output.AlgorithmOutput] = options.Algorithm(T, U0, 'MaxIter', options.MaxIter, ...
        'TolX', options.Tol, 'TolFun', options.Tol^2, 'Complex', options.Complex);
end

% r = 1; n = 1;
% U{1,r}{1,n}./U{1,r}{1,n}(1,:)

output.Options = options;

return


