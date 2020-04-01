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
eigvals = fliplr(flip(eigvals));
%% Task 2. Verify recovery
u_0z_reco = 0;
for ii=1:length(eigvals)
    u_0z_reco = u_0z_reco + eigvals(ii,ii)*sum(eigfuns{ii})*eigfuns{ii};
end
%% Task 3
lambda = diag(eigvals);
diff_eigenfunc = diff(eigfuns);
innervec = zeros(1,50);
psi_n    = innervec;
psi_dn   = innervec;
for ii=1:length(eigvals)
    innervec(ii) = sum(eigfuns{ii});
    psi_n(ii)    = eigfuns{ii}(-2*pi);
    psi_dn(ii)   = diff_eigenfunc{ii}(-2*pi);
end 
%% Task 4
fv0 = reshape([psi_dn; psi_n], [], 1);
utz = 0;
for j = 1:50
    utz = utz + lambda(j)*innervec(j)*psi_n(j);
end
utz = u_0z_reco(-2*pi);
ODEP = zeros(100);
for ii=1:50
    ODEP(2*ii,2*ii-1) = 1;
    ODEP(2*ii-1,2*ii) = lambda(ii)-innervec(ii);
end
[z2,y2] = ode113( @(z,y) myode(z,y,ODEP), 2*pi, fv0);
plot(y2(:,99))
%% Task 5 KdV numerically (as in chebfun.org)
tmax = 0.1;
S = spinop(d,[0 tmax]);
S.lin = @(u) - diff(u0 ,3);
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
for ii=1:size(lambda,1)
    Almat(ii,2*ii) = -4*lambda(ii)^2;
    Almat(2*ii,ii) = -4*eigvals(ii);
end
T = 0.1;
aa = expm(Almat)*fv0;

%% Task 7
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




%%
 