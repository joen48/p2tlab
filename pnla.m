function root = pnla(polysys, varargin)
%PNLA Wrapper for sparf.
%   root = pnla(polysys, options) solves polysys, which is a cell
%   containing the coefficients and monomial exponents of a set of
%   polynomial equations, in a similar way as PNLA_MATLAB_OCTAVE/sparf.m.

p = inputParser;
p.addOptional('degree', dbound(polysys)+1);
p.addOptional('Nsol', []);          % in a LS sense
p.parse(varargin{:});

options = cell2struct(struct2cell(p.Results), fieldnames(p.Results));

% number of variables
n = size(polysys{1,2}, 2);

if ~isempty(options.Nsol)
    m = options.Nsol;
else
    m = bezout(polysys);
end

if ~isempty(options.Nsol)
    d = options.degree;
    [M, ~] = getM(polysys, d);
    [~, ~, V] = svd(M);
    K = V(:,end-m+1:end);
else
    d = getD0(polysys);
    M = getM(polysys, d, 1);
    [Q, R, ~] = qr(M');
    r = nnz(diag(R));
    K = Q(:,r+1:end);
    while d < options.degree
        d = d + 1;
        K = updateN(K, getMex(polysys, d, d-1, 1), 1);
    end
    K = full(K);
end

db = nchoosek(n+d-1,n);
B = K(1:db,:);
A = zeros(size(B));
% better to shift with random linear combination of components in
% order to tackle multiplicities of a component
coef = randn(1,n);  % take random components
for j = 1:n
    for i = 1:db            
        A(i,:) = A(i,:) + coef(j)*K(feti([zeros(1,j-1) 1 zeros(1,n-j)] + fite(i,n)), :);
    end
end

[V, ~] = eig(pinv(B)*A); 
X = K(1:n+1,:)*V;
X = X*diag(1./X(1,:));
root = X(2:end,:)';

end
