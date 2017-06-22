function [X, output] = poly_sd(polysys, varargin)
%POLY_SD Multivariate polynomial root-finding by a SD.
%   [X, output] = poly_sd(polysys, varargin) solves polysys, which is a cell
%   containing the coefficients and monomial exponents of the set of polynomial equations.
%   X(:,k) is the root x^{(k)}.
%
%   The structure output contains the output of the simultaneous
%   diagonalization algorithm and additional information on the options used. 
%
%   Options are:
%      options.degree               If present, the degree of the Macaulay
%                                   matrix that is constructed to solve the
%                                   set of polynomial equations.
%      options.Nsol                 If present, the number of solutions,
%                                   e.g., when solving an over-constrained
%                                   set. This options will override 
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
%                                   sufficiently high. 
%                                   By default, 'auto' is used.
%      options.Complex              If true, compute a complex solution.
%      = [true|false]
%      options.Algorithm            The main algorithm for computing the CPD.
%      = [@cpd|@cpd_gevd]           If there are roots at infinity, @cpd is preferred. 
%                                   For speed, @cpd_gevd is preferred. 
%                                   By default, @cpd_gevd is used.
%
%   See also:
%       aln, getM, updateN, getRowSelect0, cpd3_sd, cpd_gevd, cp_extQZ, rescale
%
%   References:
%   [1] L. De Lathauwer. A link between the canonical decomposition in multilinear
%       algebra and simultaneous matrix diagonalization. SIAM Journal on Matrix
%       Analysis and Applications, 28(3):642--666, 2006.
%   [2] M. Sørensen and L. De Lathauwer. Multidimensional Harmonic Retrieval via
%       Coupled Canonical Polyadic Decomposition --- Part II: Algorithm and Multirate
%       Sampling. IEEE Transactions on Signal Processing, 65(2), 2017.

p = inputParser;
p.addOptional('degree', []);
p.addOptional('Nsol', []);
p.addOptional('sparse', false);
p.addOptional('RecursiveOrthogonalization', false);
p.addOptional('Complex', true);
p.addOptional('Algorithm', @cpd_gevd);
p.addOptional('Tol', 1e-6);         % should actually depend on the accuracy of the polynomial coefficients
p.addOptional('MaxIter', 500);      % 200, 500, 10000
p.parse(varargin{:});

options = cell2struct(struct2cell(p.Results), fieldnames(p.Results));

% number of equations
s = size(polysys,1);
% number of variables
n = size(polysys{1,2},2);
if (s < n)
    error('poly_sd:solution', 'The solution set is not 0-dimensional.')
end
% this vector will contain the degree of each equation
d_i = zeros(1,s);
for i = 1:s
    d_i(i) = max(sum(polysys{i,2},2)); 
end
if ~isempty(options.Nsol)
    m = options.Nsol;
    options.RecursiveOrthogonalization = false;
else
    % Bézout number
    m = prod(d_i);
    % Do not update options.Nsol for output yet. 
    % If options.Nsol is not enforced, the size of null(M(d)) will be determined automatically.
end
% % The degree of regularity is revealed at a degree \le \sum d_i + 1.
% [dp, dm] = aln(polysys, dbound + n + 1);
% d = max([dp, dm])-n;
% upper bound for the degree of regularity
dbound = sum(d_i)-n;

% if ischar(options.RecursiveOrthogonalization) && strcmpi(options.RecursiveOrthogonalization, 'auto')
%     % heuristic: 
%     % The number of rows of M(d) grows faster than the number of columns. 
%     % If the upper bound of dbound x SVD of size nrows = s*(dbound-min(d_i))
%     % is less than the upper bound of 1 x SVD of size nrows =
%     % s*nchoosek(dbound-min(d_i)+n,n), then update the null space
%     % recursively.
%     if ( (dbound-max(d_i)+1)*(s*nchoosek(n-1+dbound-min(d_i),dbound-min(d_i)))^2 ...
%             < (s*nchoosek(dbound-min(d_i)+n,n))^2 )
%         options.RecursiveOrthogonalization = true;  
%     else
%         options.RecursiveOrthogonalization = false;
%     end
% end
% if ~islogical(options.RecursiveOrthogonalization)
%     warning('poly_sd:RecursiveOrthogonalization', 'Invalid RecursiveOrthogonalization option. Will use false instead.')
%     options.RecursiveOrthogonalization = false;
% end
% options.RecursiveOrthogonalization = true;    % Krylov?

if options.sparse
    options.RecursiveOrthogonalization = false;
end

output.Time.M = cputime;
% Get numerical basis for the null space of M(d)
if options.RecursiveOrthogonalization
%     output.Time.updateN = cputime;
    d0 = max(d_i);
    K = null(getM(polysys, d0));
    if ~isempty(options.degree)
        stop = @(Z,d)(d==options.degree);
    else
        stop = @(Z,d)(size(Z,2)>=m);
    end
    d = d0;
    while ~stop(K, d)
        d = d + 1;
        [K, tol] = updateN(K, getMex(polysys, d, d-1));     % ... PART II
%         if (tol > options.Tol)    % is this the right tolerance to
%         compare with?
%             warning('poly_sd:Tol', ['The specified tolerance options.Tol = ' num2str(options.Tol) ' could not be reached.']); 
%         end
    end
%     output.Time.updateN = cputime-output.Time.updateN;
%     disp(['CPU Time for updateN = ' num2str(output.Time.updateN) ' s.'])
else
    if ~isempty(options.degree)
        d = options.degree;
    else
        d = dbound;
    end
%     output.Time.getM = cputime;
    % Get full Macaulay matrix
    [M, nnz] = getM(polysys, d);
%     output.Time.getM = cputime-output.Time.getM;
%     disp(['CPU Time for getM = ' num2str(output.Time.getM) ' s.'])
%     disp(['The relative number of nonzero elements in M(' num2str(d) ') = ' num2str(nnz/numel(M)) '.']);
%     output.Time.null = cputime;
    if ~isempty(options.Nsol)
        if options.sparse
            M = sparse(M);
            [V, ~] = spqr(M');
        else
            [~, S, V] = svd(M);
%             semilogy(diag(S), '.-', 'MarkerSize', 10)
        end
        K = V(:,end-m+1:end);
    else
        if options.sparse
            error('poly_sd:sparse', 'Cannot determine the nullity if options.sparse = true.')
        end
        K = null(M);    % Determine the nullity automatically.
    end
%     output.Time.null = cputime-output.Time.null;
%     disp(['CPU Time for null = ' num2str(output.Time.null) ' s.'])
%     if (size(K,2) < m)
%         error('poly_sd:nullity', 'The nullity of the Macaulay matrix has not converged.'); 
%     end
end
nu = size(K,2);
if (nu == 0)
    error('poly_sd:nullity', 'The nullity of the Macaulay matrix is empty.')
end
% disp(['d = ' num2str(d)])
output.Time.M = cputime-output.Time.M;

% Express shift relations
[S, T, Z] = deal(cell(1,n+1));
[S{1}, S{2}] = getRowSelect0(d, n, 1, 0);
for j = 2:n
    [~, S{j+1}] = getRowSelect0(d, n, j, 0);
end
% Complementary projectors
for j = 1:size(S,2)
   T{j} = setdiff(1:size(K,1), S{j});
   Z{j} = K(T{j},:);
end
[Y, Z] = deal(cell(1,n+1));
for j = 0:n
    Y{j+1} = K(S{j+1},:);
    Z{j+1} = K(T{j+1},:);
end

% Y = cat(1, Y{:});
E = full(cat(1, Y{:}));
% [E, ~, ~] = svd(E, 'econ');
% E = E(:,1:nu);
% for j = 0:n
%     Y{j+1} = E(j*J+1:(j+1)*J,:);
% end

% Two-dimensional Harmonic Retrieval BY SIMULTANEOUS MATRIX DIAGONALIZATION
% Compute the Gramian P'*P for the rank-one detecting device P. 
% The Gramian is computed as a fourth-order tensor in which element with indices
% (r,s,t,u) is equal to <Er,Et>*<Es,Eu> + <Er,Eu>*<Es,Et> - <Er'*Et,Es'*Eu>
% - <Er'*Eu,Es'*Et>, where En denotes the matricized n-th column of E.

output.Time.T = cputime;
% output.Time.sd = cputime;
YY = cell(n+1,n+1);
E12 = zeros(nu);
for j = 0:n
%     YY{j+1,j+1} = Y{j+1}'*Y{j+1};
    YY{j+1,j+1} = eye(nu) - Z{j+1}'*Z{j+1};
    E12 = E12 + YY{j+1,j+1};
    for l = j+1:n
        YY{j+1,l+1} = Y{j+1}'*Y{l+1};
    end
end
% Compute all products P12(r,s,t,u) = <Er,Et>*<Es,Eu>
P12 = reshape(kron(E12, E12), [nu nu nu nu]);
% Compute all products P34(r,s,t,u) = conj(<Er'*Et,Es'*Eu>).
P34 = zeros(nu^2, (n+1)^2);
X340 = cat(2, YY{1,:});  
P34(:,1:(n+1)) = reshape(X340, nu^2, n+1);
for j = 1:n
    X = cellfun(@ctranspose, YY(1:j,j+1), 'UniformOutput', false);
    X = [ cat(2, X{:}) cat(2, YY{j+1,j+1:n+1}) ];
    P34(:,j*(n+1)+(1:(n+1))) = reshape(X, nu^2, n+1);
end

P34 = reshape(conj(P34)*(P34.'), [nu nu nu nu]);
% Compose the matrix P'*P,
PP = reshape(permute(P12-P34, [1 4 2 3]), [nu^2 nu^2]) ...
   + reshape(permute(P12-P34, [1 4 3 2]), [nu^2 nu^2]);
% without redundant columns/rows.
lt = tril(true(nu,nu));
PP = PP(lt(:),lt(:));
% Compute the kernel of the adjusted Gramian and build the tensor M,
% consisting of symmetric matrices Mk = W Lk W^T.
[~, ~, V] = svd(PP, 'econ');
% [V, ~] = qr(PP');
V = V(:,end:-1:end-nu+1);
rr = zeros(1,nu*(nu+1)/2);
rr([1 1+cumsum(nu:-1:2)]) = 1;
V(~rr,:) = 0.5*V(~rr,:);
M = zeros(nu,nu,nu);
Mlt = zeros(nu,nu);
for r = 1:nu
    Mlt(lt) = V(:,r);
    M(:,:,r) = Mlt + tril(Mlt,-1).';
end

if ~islogical(options.Complex)
    warning('poly_sd:Complex', 'Invalid Complex option. Will use true instead.')
    options.Complex = true;
end 
options.IsReal = ~options.Complex;

% Simultaneously diagonalize the frontal slices of M.
options.AlgorithmOptions = struct;
options.AlgorithmOptions.MaxIter = options.MaxIter;
options.AlgorithmOptions.TolX = options.Tol;
options.AlgorithmOptions.TolFun = options.Tol^2;
options.AlgorithmOptions.Complex = options.Complex;
options.AlgorithmOptions.IsReal = options.IsReal;
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if ~xsfunc(options.Algorithm)
    error('poly_sd:Algorithm', 'Invalid Algorithm option.')
end
% [output.iter, W, B, C] = cp_extQZ(M, [options.MaxIter, options.Tol^2, 1]);
if strcmpi(func2str(options.Algorithm), 'cpd_gevd') 
    [W, output.AlgorithmOutput] = cpd_gevd(M, m, options.AlgorithmOptions);    % cannot handle roots at infinity properly
    W = W{1};
    % Here, we are assuming $x_0 = 1$:
    den = diag(W'*X340(:,1:nu)*W);
    X = [ones(1,m); zeros(n,m)];
    for j = 1:n 
       X(j+1,:) = diag(W'*X340(:,j*nu+(1:nu))*W)./den;
    end
elseif strcmpi(func2str(options.Algorithm), 'cpd') 
    [W, output.AlgorithmOutput] = cpd(M, m, options.AlgorithmOptions, 'Initialization', @cpd_rnd);
    W = W{1};
    % This works for the non-generic case as well.
    krAB = conj(E*W);
    X = zeros(n+1,m);
    for r = 1:m
        [~, ~, v] = svd(reshape(krAB(:,r), [], n+1), 'econ');
        X(:,r) = conj(v(:,1));
    end
    X = rescale({X, zeros(1, size(X,2)), zeros(1, size(X,2))});     % A BETTER VERSION OF RESCALE?
    X = X{1};
else
    error('poly_sd:Algorithm', 'Invalid Algorithm option.')
end
% output.Time.sd = cputime-output.Time.sd;
output.Time.T = cputime-output.Time.T;

output.Options = options;

end

