% SC42135 SPECTRAL ANALYSIS OF NONLINEAR/INFINITE-DIMENSIONAL SYSTEMS
% Practice session 3 Nonlinear Fourier transforms
clear all; close all; clc;
%% Task 1. Eigenvalues of the Schröedinger operator
d = [-2*pi,2*pi];
A = chebop(d); 
A.bc = 'periodic';
x = chebfun('x', d);
u0 = 4*exp(-x^2);
A.op = @(u) diff(u,2)+u0*u;
[eigfuns,eigvals] = eigs(A,50);
max(eigvals(:))
%% Task 2. Verify recovery
u_0z = eigvals(1,1)*eigfuns{1}^2;
for ii=2:length(eigvals)
    u_0z = eigvals(ii,ii)*eigfuns{ii}^2;
end
%% Task 3
lambda = diag(eigvals);
diff_eigenfunc = diff(eigfuns);
innervec = sum(eigfuns{1}, -2*pi, 2*pi);
for ii=2:length(eigvals)
    innervec = sum(eigfuns{ii}, -2*pi, 2*pi);
end 
%% Task 4
fv0 = reshape([diff_eigenfunc; eigfuns], [], 1);


%%[TOUT,YOUT] = ode113(ODEFUN,TSPAN,Y0);
%% Task 5 KdV numerically (as in chebfun.org)
tmax = 0.1;
S = spinop(d,[0 tmax]);
S.lin = @(u) - diff(u,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.nonlin = @(u) -3*diff(u.^2); % spin cannot parse "u.*diff(u)"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.init = u0;
N = 800;   % numer of grid points
dt = 5e-6; % time-step
u_t01 = spin(S,N,dt,'plot','off');
%% Task 6 
A_t01 = chebop(d); 
A_t01.bc = 'periodic';
A_t01.op = @(u) diff(u,2)+u_t01*u;
[eigfunsComp,eigvalsCom] = eigs(A_t01,50);
max(eigvalsCom(:))

Almat = zeros(N*2);
for ii=1:size(eigvals,1)
    Almat(1,2) = -4*eigvals(1,1)^2;
    Almat(2,1) = -4*eigvals(1,1);
end

expm()

%% Task 7
%% Plots
figure;

% task 5
plot(S.init), hold on, plot(u), hold off
text(4.4,1300,'t = 0'), text(13.5,1300,'t = 0.0156')