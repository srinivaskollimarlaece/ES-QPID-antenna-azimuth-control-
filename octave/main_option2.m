% MAIN_OPTION2  Complete ES-QPID study: all experiments and figures of the article.
%
%   Runs every experiment (about 10 to 15 minutes in GNU Octave 8), prints all tables,
%   writes results_option2.txt / results_option2.mat and saves the figures.
%   Requires the Octave control package (pkg load control). For MATLAB use the folder ../matlab.

clc; clear; close all;
pkg load control
tic;
R = sqpid_experiments(2, 200);      % 200 Monte Carlo runs per delay
sqpid_figures(2);
printf('Finished in %.0f s\n', toc);
