function C = compoundk(A, k)
%COMPOUNDK The k-th compound matrix of A
%   The k-th compound matrix of A is denoted by Ck.
%   It is the matrix containing the determinants of all k xk submatrices of A,
%   arranged with the submatrix index sets in lexicographic order.

r = nchoosek(1:size(A,1), k);
c = nchoosek(1:size(A,2), k);
C = zeros(size(r,1), size(c,1));

for i = 1:size(C,1)
    for j = 1:size(C,2)
        C(i,j) = det(A(r(i,:),c(j,:)));
    end
end

