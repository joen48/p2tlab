function U = rescale(U, Tol)
%RESCALE Rescale poly_sd solution.
%   U = rescale(U, Tol) removes scaling ambiguity from the CPD
%   solution columns, such that U{1}(1,:) == 1 if U{1}(1,:) >= Tol originally,
%   or U{1}(m(k),k) == 1 where m(k) is the position of the maximal element
%   in the kth column originally otherwise.
%
% See also:
%   poly_cpd, poly_sd

if (nargin < 2), Tol = 1e-6; end;

N = size(U, 2);
for r = 1:size(U{1}, 2)
    idx = 1; 
    if (abs(U{1}(idx,r)) < Tol*norm(U{1}(:,r)))
       [~, idx] = max(abs(U{1}(:,r)));
    end 
    lambda = U{1}(idx,r);
    U{1}(:,r) = U{1}(:,r)/lambda;   
    U{N}(:,r) = lambda*U{N}(:,r);    
end

end

%RESCALE Rescale MHR factors
%   U = rescale(U) rescales the MHR factor matrices U{1}, ..., U{N-1}
%   such that U{1:N-1}(1,:) = 1
% 
% N = size(U, 2);
% for r = 1:size(U{1},2)
%     lambda = zeros(1, N-1);
%     for n = 1:length(lambda)   
%         lambda(n) = U{n}(1,r);
%         U{n}(:,r) = U{n}(:,r)/lambda(n);   
%     end 
%     U{N}(:,r) = prod(lambda)*U{N}(:,r);    
% end
% 
% end
