clear, clc, close all

results = runtests('cpd_test.m');
rt = table(results);
sortrows(rt, 'Duration')

% results = runtests('sd_test.m');
% rt = table(results);
% sortrows(rt, 'Duration')

% results = runtests('merge_test.m');
% rt = table(results);
% sortrows(rt, 'Duration')