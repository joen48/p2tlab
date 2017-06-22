function tests = merge_test
tests = functiontests(localfunctions);
end

% function setupOnce(testCase)
% set_phcpath([pwd filesep 'phc']);
% end

function setup(testCase)
rng(0, 'v5uniform');
end

function test114(testCase)
% Example 1.1.4
% Back to the Roots (2.1)
% system of 2 equations in 2 variables with 4 affine solutions
polysys = {};
polysys{1,1} = [-4 5 -3 -1 2 1]';
polysys{1,2} = [0 1 0 2 1 0; 0 0 1 0 1 2]';
polysys{2,1} = [-1 0 0 1 2 1]';
polysys{2,2} = [0 1 0 2 1 0; 0 0 1 0 1 2]';
X = [1 1 1 1; 0 1 3 4; -1 0 -2 -5];
options = struct;
options.Tol = 1e-6;
[Xest, ~] = mergesolve(polysys, options);
relerr = cpderr(X, Xest);
verifyLessThan(testCase, relerr, options.Tol);
end