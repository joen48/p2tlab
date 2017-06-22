function [X, output] = poly_cpd(polysys, varargin)
%POLY_CPD Multivariate polynomial root-finding by a CPD.
%   [X, output] = poly_cpd(polysys, varargin) solves polysys, which is a cell
%   containing the coefficients and monomial exponents of the set of polynomial equations. 
%   X(:,k) is the root x^{(k)}.
%
%   The structure output contains the output of the algorithm used for
%   computing the CPD and additional information on the options used. 
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
%      options.Algorithm            The main algorithm for computing the CPD.
%      = [@cpd||@cpd_gevd|@cpd_nls] By default, @cpd is used.
%
%   See also:
%       aln, getM, updateN, getRowSelect0, cpd, rescale

p = inputParser;
p.addOptional('degree', []);
p.addOptional('Nsol', []);
p.addOptional('sparse', false);
p.addOptional('RecursiveOrthogonalization', false);
p.addOptional('Complex', true);
p.addOptional('Compression', 'svd');
p.addOptional('Algorithm', @cpd);
p.addOptional('Tol', 1e-6);         % should actually depend on the accuracy of the polynomial coefficients
p.addOptional('MaxIter', 500);      % 200, 500, 10000
p.parse(varargin{:});

options = cell2struct(struct2cell(p.Results), fieldnames(p.Results));

% number of equations
s = size(polysys, 1);
% number of variables
n = size(polysys{1,2}, 2);
if (s < n)
    error('poly_ll1:solution', 'The solution set is not 0-dimensional.')
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
    % % The degree of regularity is revealed at a degree \le \sum d_i + 1.
    % [dp, dm] = aln(polysys, dbound + n + 1);
    % d = max([dp, dm])-n;
% upper bound for the degree of regularity
dbound = sum(d_i)-n;

% @cpd_gevd fails if d < d*+1 
if strcmpi(func2str(options.Algorithm), 'cpd_gevd')
    options.degree = dbound+1;  % Calculate a numerical basis for the null space of M(d*+1)
end
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
if strcmpi(options.Compression, 'none') && strcmpi(func2str(options.Algorithm), 'cpd_gevd')
    warning('poly_cpd:Compression', 'Compression is needed with option.Algorithm = @cpd_gevd. Will use ''svd''.')
    options.Compression = 'svd';
end
%     if strcmpi(options.Compression, 'mon')  % not correct
%     %     [~, b, ~] = candecomp(polysys, d);
%         b = zeros(m, 1);
%         numel_b = 0;
%         for k = 1:size(K,1)
%             if (numel_b == m), break; end;
%             if (rank(K([b(1:numel_b); k],:)) > numel_b)
%                 numel_b = numel_b+1;
%                 b(numel_b) = k;
%             end
%         end
%         b_low = b(b <= S{1}(end)); % standard monomials up to degree d-1
%         J = numel(b_low);
%         for j = 0:n
%             Y{j+1} = K(S{j+1}(b_low),:);
%         end

T = polysys2tens(polysys, 1, struct('degree', options.degree, ...
    'm', m, 'sparse', options.sparse, 'RecursiveOrthogonalization', options.RecursiveOrthogonalization, ...
    'd0', max(d_i), 'Compression', options.Compression));

if (size(T,3) == 1)    % only 1 solution, so it is readily available
    U = rescale({T(:,:,1), 1, 1}, options.Tol);
    X = U{1};
    output.Options = options;
    return
end

% CPD
options.AlgorithmOptions = struct;
if ~islogical(options.Complex)
    warning('poly_cpd:Complex', 'Invalid Complex option. Will use true instead.')
    options.Complex = true;
end 
options.IsReal = ~options.Complex;
options.AlgorithmOptions.MaxIter = options.MaxIter;
options.AlgorithmOptions.TolX = options.Tol;
options.AlgorithmOptions.TolFun = options.Tol^2;
options.AlgorithmOptions.Complex = options.Complex;
options.AlgorithmOptions.IsReal = options.IsReal;
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if ~xsfunc(options.Algorithm)
    error('poly_cpd:Algorithm', 'Invalid Algorithm option.')
end
% output.Time.Algorithm = cputime;    
if strcmpi(func2str(options.Algorithm), 'cpd') || strcmpi(func2str(options.Algorithm), 'cpd_gevd')
    [U, output.AlgorithmOutput] = options.Algorithm(T, m, options.AlgorithmOptions);
else
    options.AlgorithmOptions.LargeScale = (sum(size(T))*m > 1e2);
    if options.Complex
        U0 = cpd_rnd(size(T), m, 'Real', @randn, 'Imag', @randn);
    else
        U0 = cpd_rnd(size(T), m, 'Real', @randn);
    end
    [U, output.AlgorithmOutput] = options.Algorithm(T, U0, options.AlgorithmOptions);
end
% output.Time.Algorithm = cputime-output.Time.Algorithm;

% Visualize convergence
% figure
% semilogy(output.AlgorithmOutput.relfval), hold on
% semilogy(1:output.iterations+1, options.Tol^2*ones(1,output.iterations+1), 'k--')
% ylabel('Objective function'), xlabel('Iteration')
% title('Convergence plot')

% figure(10), imagesc(log10(abs(U{2}))), colorbar

% Multiplicities?
% [3] A. Stegeman. A three-way Jordan canonical form as limit of low-rank tensor
%     approximations. SIAM Journal on Matrix Analysis and Applications, 34(2):524?
%     -650, 2013.
%     relerr = Inf;
%     U = {};
%     for itrial = 1:10   % For each array, we run CP ALS 10 times with random starting values
%         [Utrial, out] = cpd(T, m, 'Initialization', @cpd_rnd, 'Algorithm', @cpd_als, 'MaxIter', 10000, 'TolFun', 1e-9);     % We use convergence criterion 1e-9 in CP ALS.
%         if (out.Algorithm.relerr < relerr)
%             relerr = out.Algorithm.relerr;
%             U = Utrial;
%         end
%     end
%     Ut = kr(U);
%     UU = Ut'*Ut;
%     d = sqrt(diag(UU));
%     UU = UU./(d*d');                % equal to the cosine of the angle between the vectorized rank-1 terms
%     lt = tril(abs(UU) > 0.95);      
%     g = {};
%     for s = 1:m
%         ng = setdiff(find(lt(:,s)'), [g{:}]);   % mind transitivity
%         if ~isempty(ng)
%             g{end+1} = ng;
%         end
%     end

% if (cellfun(@numel, g) == 1)    % no multiplicities > 1
    U = rescale(U, options.Tol);
%     norm(diag(cpd_crb(U, 1, 'Method', 'full')));
    X = U{1};
    output.Options = options;
% else
%     warning('poly_cpd:multiplicities', 'Some roots have multiplicity greater than 1.')
% end
% disp(['CPU Time for cpd = ' num2str(toc) ' s.']) 
return

% Theorem 3.2.1
C2 = kr(compoundk(U{1}, 2), compoundk(U{2}, 2));
if (rank(U{3}) < m) || (rank(C2) < size(C2, 2)) 
    warning('poly_cpd:uniqueness', 'The solution might not be unique.')
end
return