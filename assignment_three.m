% SC42135 SPECTRAL ANALYSIS OF NONLINEAR/INFINITE-DIMENSIONAL SYSTEMS
% Practice session 3 Nonlinear Fourier transforms
clear all; close all; clc;
%% Task 1. Eigenvalues of the Schröedinger operator

d = [-2*pi,2*pi];
A = chebop(d); 
A.bc = 'periodic';
x = chebfun('x', d);
u0 = 4*exp(-x^2);
A.op = @(u) +diff(u,2)+u0*u;
[vv,D] = eigs(A,50);
 max(D(:))
