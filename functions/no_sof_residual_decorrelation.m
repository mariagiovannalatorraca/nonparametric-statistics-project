% DISCLAIMER:
% The matrix sigma_res is symmetric and positive definite, allowing us
% to compute its inverse through LU factorization and then utilize
% Cholesky decomposition to calculate its square root.
clear all
clc

% Carico la matrice sigma_res
res = csvread('data\no_sof_residuals.csv')';
sigma_res = csvread('data\no_sof_sigma.csv');

% Calcolo i residui decorrelati e li salvo in un csv
chol_sigma_res = chol(sigma_res);
csvwrite('data\no_sof_squared_root_sigma.csv', chol_sigma_res);
res_decorr = chol_sigma_res\res;
csvwrite('data\no_sof_decorrelated_residuals.csv', res_decorr);


