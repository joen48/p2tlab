function tests = sd_test
tests = functiontests(localfunctions);
end

% function setupOnce(testCase)
% set_phcpath([pwd filesep 'phc']);
% end

function setup(testCase)
rng(0, 'v5uniform');
end

function test111(testCase)
% Example 1.1.1
% system of 2 equations in 2 variables with 2 affine solutions with
% multiplicity 1 and 1 affine solution with multiplicity 2 
polysys = {};
polysys{1,1} = [1 -2]';
polysys{1,2} = [1 0; 1 1]';
polysys{2,1} = [2 -1]';
polysys{2,2} = [0 2; 2 0]';
X = [repmat([1; 0; 0], [1 2]), [1 1; 2 2; sqrt(2) -sqrt(2)]];
options = struct;
options.Tol = 1e-6;
[Xest, ~] = poly_sd(polysys, options);
relerr = cpderr(X, Xest);
verifyLessThan(testCase, relerr, options.Tol);
end

function unitest(testCase)
% Example 1.1.2
% s=1 equation in n=1 variable
polysys = {};
polysys{1,1} = [6 -5 1]';
polysys{1,2} = [2 1 0]';
X = [1 1; 1/2 1/3];
options = struct;
options.Tol = 1e-6; options.RecursiveOrthogonalization = true;
[Xest, ~] = poly_sd(polysys, options);
relerr = cpderr(X, Xest);
verifyLessThan(testCase, relerr, options.Tol);
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
[Xest, ~] = poly_sd(polysys, options);
relerr = cpderr(X, Xest);
verifyLessThan(testCase, relerr, options.Tol);
end

function test214(testCase)
% Example 2.1.4
% system of 2 equations in 2 variables with 1 affine solution and 1 root
% at infinity
polysys = {};
polysys{1,1} = [1 -2]';
polysys{1,2} = [1 0; 1 1]';
polysys{2,1} = [1 -3]';
polysys{2,2} = [0 0; 1 0]';
X = [1 0; 2 1; 3 0];
options = struct;
options.Tol = 1e-8; options.Algorithm = @cpd;
[Xest, ~] = poly_sd(polysys, options);
relerr = cpderr(X, Xest);
verifyLessThan(testCase, relerr, options.Tol);
end

function testUniquee(testCase)
% Example -
% system of 2 equations in 2 variables with 2 affine solutions and 2 roots
% at infinity
polysys = {};
polysys{1,1} = [1 -1]';
polysys{1,2} = [1 1; 1 0]';
polysys{2,1} = [1 -1]';
polysys{2,2} = [1 0; 1 1]';
X = [1 1 0 0; 0 1 1 0; 0 1 0 1];
options = struct;
options.Tol = 1e-6; options.Algorithm = @cpd;
[Xest, ~] = poly_sd(polysys, options);
relerr = cpderr(X, Xest);
verifyLessThan(testCase, relerr, options.Tol);
end

function test331(testCase)
% Example 3.3.1
% 145 Example 7
% system of 2 equations in 2 variables with 1 affine solution with
% multiplicity 4
polysys = {};
polysys{1,1} = [1 -4 4]';
polysys{1,2} = [0 0 0; 2 1 0]';
polysys{2,1} = [1 2 -2 1 -2 1]';
polysys{2,2} = [0 1 0 2 1 0; 0 0 1 0 1 2]';
X = repmat([1; 1; 2], [1 4]);
options = struct;
options.Tol = 1e-6; options.RecursiveOrthogonalization = true; options.Algorithm = @cpd_gevd;
[Xest, ~] = poly_sd(polysys, options);
relerr = cpderr(X, Xest);
verifyLessThan(testCase, relerr, options.Tol)
end

function test3(testCase)
% Back to the Roots (2.3)
% system of 3 equations in 3 variables with 18 affine solutions
polysys = {};
polysys{1,1} = [1 -1 1]';
polysys{1,2} = [2 1 0; 0 1 0; 0 0 1]';
polysys{2,1} = [1 -2 -3]';
polysys{2,2} = [0 1 1; 3 2 1; 0 0 0]';
polysys{3,1} = [1 -1 -2]';
polysys{3,2} = [0 1 0; 0 1 0; 3 1 0]';
load('23.mat')
options = struct;
options.Tol = 1e-6; options.RecursiveOrthogonalization = true;
[Xest, ~] = poly_sd(polysys, options);
relerr = cpderr(X, Xest);
verifyLessThan(testCase, relerr, options.Tol)
end
