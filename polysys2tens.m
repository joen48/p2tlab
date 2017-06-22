function T = polysys2tens(polysys, L, varargin)
%POLYSYS2TENS Convert a set of polynomial equations into a third-order tensor.
%   T = polysys2tens(polysys, options) uses spatial smoothing in every mode
%   [1; x^(j), x^(j)2, ..., x^(j)L] to convert polysys into a third-order tensor.
%
%   Options are:
%      options.degree               The degree of the Macaulay matrix. By
%                                   default, the degree of regularity + L.
%      options.m                    The size of null(M(d)). By default, the
%                                   Bézout number.
%      options.sparse               If true, compute the numerical basis for 
%                                   the null space of the Macaulay matrix 
%                                   using a sparse QR algorithm.
%      options.RecursiveOrthogonalization If true, update the numerical
%      = [true|false]               basis for the null space of the Macaulay
%                                   matrix using the recursive
%                                   orthogonalization theorem. By default,
%                                   false.
%      options.d0                   If options.RecursiveOrthogonalization is true,
%                                   the degree to start constructing the
%                                   Macaulay matrix.
%      options.Compression          Compress the tensorized null space of
%      = ['svd'|'none']             the Macaulay matrix by means of an
%                                   (ML)SVD. By default, 'svd' is used.
%
%   See also:
%       poly_cpd, poly_btd

p = inputParser;
p.addOptional('degree', []);
p.addOptional('m', bezout(polysys));
p.addOptional('sparse', false);
p.addOptional('RecursiveOrthogonalization', false);
p.addOptional('d0', []);
p.addOptional('Compression', 'svd');
p.parse(varargin{:});

options = cell2struct(struct2cell(p.Results), fieldnames(p.Results));

n = size(polysys{1,2},2);
m = options.m;

if isempty(options.degree)
    options.degree = dbound(polysys) + L;
else
    warning('polysys2tens:degree', 'You specified the degree of the Macaulay matrix. Make sure that it is consistent with L.')
end

% Get numerical basis for the null space of M(d)
if options.RecursiveOrthogonalization
    if ~isfield(options, 'd0')
        error('polysys2tens:d0', 'If options.RecursiveOrthogonalization, options.d0 should be set.')
    end
    d0 = options.d0;
    K = null(getM(polysys, d0));
    d = d0;
    stop = @(Z,d)(d>=options.degree);
    while ~stop(K, d)
        d = d + 1;
        [K, ~] = updateN(K, getMex(polysys, d, d-1));
%         if (tol > options.Tol)    % is this the right tolerance to
%         compare with?
%             warning('poly_cpd:Tol', ['The specified tolerance options.Tol = ' num2str(options.Tol) ' could not be reached.']); 
%         end
    end
else
    d = options.degree;
    % Get full Macaulay matrix
    M = getM(polysys, d);
    if options.sparse
        M = sparse(M);
        [V, ~] = spqr(M');
    else
        [~, ~, V] = svd(M);
    end
    K = V(:,end-m+1:end);
end

% Express shift relations
S = cell(1,1+L*n);
[S{1}, ~] = getRowSelect0(d-L+1, n, 1, 0);
for j = 1:n
    [~, Sj] = getRowSelect0(d, n, j, 0);    % spatially smooth once by default
    S{1+(j-1)*L+1} = Sj(S{1});
    for l = 2:L
        S{1+(j-1)*L+l} = Sj(S{1+(j-1)*L+l-1});
    end
end
Y = cell(size(S));
for i = 1:size(Y, 2)
    Y{i} = K(S{i},:);
end

if ~strcmpi(options.Compression, 'none') && (size(S{1}, 2) < m)
    warning('polysys2tens:Compression', 'No compression is needed.')     % should be mutually exclusive with the previous condition
end
if strcmpi(options.Compression, 'svd')
    Z = cat(2, Y{:});
    if options.sparse
        [Uz, Sz, ~] = svds(Z, m);   % issparse(spqr) = true
    else
        [Uz, Sz, ~] = svd(Z);
    end
    tol = max(size(Z)) * eps(max(diag(Sz))); J = min(sum(diag(Sz) > tol), m);   % avoid rank(), which calculates SVD again
    Uz = Uz(:,1:J);
    for j = 0:n
        Y{j+1} = Uz'*Y{j+1};    % equivalent to tmprod(T, Uz', 2);
    end
end

Y = cat(1, Y{:});
T = mat2tens(Y, [1+L*n size(Y,1)/(1+L*n) size(Y,2)], [2 1], 3); % size(Y,2) == m   

end

