% MAIN_OPTION2  Complete ES-QPID study: enhanced Smith predictor-based quantum-inspired PID
%   (ES-QPID) for antenna azimuth control. Runs all experiments and figures of the article.
%
%   Runs every experiment (E1 ... E7), prints all tables, writes results_option2.txt and
%   results_option2.mat, and saves every figure as PNG in the folder figures_option2.
%   Base MATLAB only (no toolbox needed). Keep sqpid_random.mat in the same folder.
%   Run time: a few minutes; most of it is the Monte Carlo study (E6).

clc; clear; close all;
tic;
R = sqpid_experiments(2, 200);      % 2 = Option 2, 200 Monte Carlo runs per delay
sqpid_figures(2);
fprintf('Finished in %.0f s\n', toc);
