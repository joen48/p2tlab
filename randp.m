function [polysys, X] = randp(d, n, varargin)
%RANDP Pseudorandom set of polynomial equations.
%   [polysys, X] = randp(d, n, varargin) returns a set of polynomial
%   equations of total degree d in n variables with normally distributed
%   pseudorandom coefficients.
%
%   Options are:
%      options.Complex              If true, the pseudorandom coefficients
%                                   are complex. By default, the
%                                   pseudorandom coefficients are real.
%      options.seed                 If present, the non-negative integer
%                                   seed to seed the random number
%                                   generator used by RANDN.
%      options.solveit              If true, also the solution X is
%                                   returned. If false, X is empty. 
%
%   See also:
%       randn, getMonBase, getM, solve_system, repset

p = inputParser;
p.addOptional('Complex', false);
p.addOptional('seed', []);
p.addOptional('solveit', false);
p.parse(varargin{:});

options = cell2struct(struct2cell(p.Results), fieldnames(p.Results));

if ~isempty(options.seed)
    rng(options.seed, 'v5uniform');
end

mon = getMonBase(d, n+1);
q = size(mon,1);
polysys = cell(n,2);
t = zeros(n*(q+1),n+1);
for i = 1:n
    polysys{i,1} = randn(size(mon,1), 1) + options.Complex*1j*randn(size(mon,1), 1);
    polysys{i,2} = mon(:,2:end);
    t((i-1)*(q+1)+1:i*(q+1)-1,:) = [polysys{i,1} polysys{i,2}];
end

X = [];
% Solve with PHClab
if options.solveit
    warning('randpoly:solveit', 'If function call returns with error, make sure PHClab is installed and set up.')
    s = solve_system(t);
    m = size(s,2);
    X = zeros(n+1, m);
    X(1,:) = ones(1,m);
    for i = 1:n
        X(i+1,:) = [s.(sprintf('x%d',i))];
    end
end