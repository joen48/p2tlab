function [Sa, Sb] = getRowSelect0(d, n, a, b)
% [Sa Sb] = getRowSelect(d,n,a,b)
% - 2017/01/04 JV

if (a < 0) || (a > n) || (b < 0) || (b > n)
    error('getRowSelect0:shift', ['Shift variables should be scalars between ' num2str(0) ' and ' num2str(n) '.'])
end

mons = getMonBase(d, n+1);
Sa = find(mons(:,b+1) > 0)';
Sb = mons(Sa,:);
Sb(:,b+1) = Sb(:,b+1) - 1;
Sb(:,a+1) = Sb(:,a+1) + 1;
Sb = feti(Sb(:,2:end));

% assert(max(Sb) <= size(mons,1));

end

