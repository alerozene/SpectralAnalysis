% SC42135 SPECTRAL ANALYSIS OF NONLINEAR/INFINITE-DIMENSIONAL SYSTEMS
% Practice session 3 Nonlinear Fourier transforms
clear all; close all; clc;

%% Task 1. Eigenvalues of the Schröedinger operator
d = [-2*pi,2*pi];
numeigs = 50;
A = chebop(d); 
A.bc = 'periodic';
x = chebfun('x', d);
u0 = 4*exp(-x^2);
A.op = @(u) diff(u,2)+u0*u;
[eigfuns,eigvals] = eigs(A,numeigs);
eigfuns = fliplr(eigfuns);
eigvals = diag(fliplr(flip(eigvals)));
%% Task 2. Verify recovery
u_0_reco = 0;
for ii=1:length(eigvals)
    u_0_reco = u_0_reco + eigvals(ii)*sum(eigfuns{ii})*eigfuns{ii};
end
%% Task 3

diff_eigfuns = diff(eigfuns);

innervec = zeros(1,50);
psi_n    = innervec;
psi_dn   = innervec;
for ii=1:length(eigvals)
    innervec(ii) = sum(eigfuns{ii});
    psi_n(ii)    = eigfuns{ii}(-2*pi);
    psi_dn(ii)   = diff_eigfuns{ii}(-2*pi);
end

%% Task 4
fv0 = reshape([psi_dn; psi_n], [], 1);
utz = 0;
for j = 1:50
    utz = utz + eigvals(j)*innervec(j)*psi_n(j);
end
utz = u_0_reco(-2*pi);
ODEP = zeros(100);
for ii=1:50
    ODEP(2*ii,2*ii-1) = 1;
    ODEP(2*ii-1,2*ii) = eigvals(ii)-innervec(ii);
end
[z2,y2] = ode113( @(z,y) P3T4(z,y,eigvals,innervec), [-2*pi 2*pi], fv0);

%% Task 5 KdV numerically (as in chebfun.org)
tmax = 0.1;
S = spinop(d,[0 tmax]);
S.lin = @(u) - diff(u ,3);
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

Almat = zeros(50*2); 
for ii=1:size(eigvals,1)
    Almat(ii,2*ii) = -4*eigvals(ii)^2;
    Almat(2*ii,ii) = -4*eigvals(ii);
end
T = 0.1;
 
= expm(T*Almat)*fv0;

%% Task 7 Compute eigvalues/functions of 3.6 and verify 3.5
A_var1 = 12;
P = 16*pi;
d_P = [-P/2,P/2];
% Repetition of task 1 changing u0 for u0_z and the domains
Atrace = chebop(d_P); 
Atrace.bc = 'periodic';
z = chebfun('z', d_P);
u0_z = A_var1*sech(z)^2;
Atrace.op = @(u) diff(u,2)+u0_z*u;
[trace_psi,trace_lambda] = eigs(Atrace,numeigs);
trace_psi = fliplr(trace_psi);
trace_lambda = diag(fliplr(flip(trace_lambda)));
trace_lambda = trace_lambda(trace_lambda>0);
trace_psi = trace_psi(1:length(trace_lambda));
% Verify 3.4
u0_zrecon = 0;
for ii=1:length(trace_lambda)
    u0_zrecon = u0_zrecon + sqrt(trace_lambda(ii))*trace_psi{ii}.^2;
end
figure
plot(u0_z)
hold on
plot(u0_zrecon)

%% Task 8 A=8 P = 16pi non longer multisoliton


%% Plots
% Task 1
figure(1);
title('Eigenfunctions of the Schrödinger operator')
hold on
for ii=1:10
    txt = ['\lambda = ',num2str(eigvals(ii,ii))];
    plot(eigfuns{ii},'DisplayName',txt)
end
hold off
legend show
% Task 2
figure(2);
title('Reconstruction u(0,z) = 4e^{-z^2}')
subplot(1,2,1)
plot(u0)
title('chebfun')
subplot(1,2,2)
plot(u_0z_reco)
title('reconstruction')


% task 5
plot(S.init), hold on, plot(u), hold off
text(4.4,1300,'t = 0'), text(13.5,1300,'t = 0.0156')