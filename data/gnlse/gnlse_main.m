% Run all other scripts for the Generalized Nonlinear Schrodinger equation
% Jared Callaham (2020)

clear all; close all; clc

%% Simulate supercontinuum generation
% Code by J.C. Travers, M.H Frosz and J.M. Dudley (2009)
% www.scgbook.info
test_Dudley;

%% Nondimensionalize results according to soliton scaling
nondimensionalize;

%% Separate terms in GNLSE
gnlse_terms;