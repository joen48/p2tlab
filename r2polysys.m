function polysys = r2polysys(r)
%R2POLYSYS Polysys with specified roots.  
polysys = {};
polysys{1,1} = poly(r)';
polysys{1,2} = (numel(r):-1:0)';
end

