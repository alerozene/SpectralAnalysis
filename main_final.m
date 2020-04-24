% SC42135 Spectral Analysis of Nonlinear/Infinite-Dimensional Systems.
% Final assignment

close all; clear all; clc
%% Part 1 NLFT
% T1&2 Compute psi(a) with chebfun and find coefficients gamma and delta
a = 2;
u1 = 1;
lambda1 = -1;

[gamma1,psi1] = eigenfunctions(a,u1,lambda1);

% Task 3 Compute|gamma| for some continuous spectrum with u=20*exp(-z^4)
tic 
a = 25;
lambda = -25:0.1:0;
zz = chebfun('z', [-a a]);
u = 20*exp(-zz^4);
[gammadelta,psivec] = eigenfunctions(a,u,lambda);
toc

%% Task 4.Stochastic dynamical model with DMD
% 3.2.1 stochastic dynamics system
lambdaslope = 0.5;
n_variance = sqrt(10);
k = 1e3;
nk = n_variance*randn([1 k]);
z = zeros([1 k]);
% eq(18)
for ii = 2:k
    z(ii) = lambdaslope*z(ii-1) +nk(ii);
end

% Domains for Z and Zp
dom1 = 1:k-1;
dom2 = 2:k;

% Data matrices
Z = z(1:k-1);
Zp = z(2:k);

% 2.2
X = Z;
Y = Zp;

% Definition 1 (lamda hat) 
A = Y*pinv(X);

% line fit
xdom = -20:0.1:20;
dmdfit = polyfit(Z,Zp,1);

%% Figures
figure(1)
subplot(211)
plot(-25:0.1:0,abs(gammadelta(1,:)),'b-','LineWidth',1.3)
subplot(212)
title('\gamma(a=25) for \lambda')
plot(-25:0.1:0,abs(gammadelta(1,:)),'b-','LineWidth',1.3)
xlabel('\lambda')
ylabel('|\gamma(\lambda)|')
grid on
axis([-10 -5 0.1 0.16])

figure(4)
subplot(121)
plot(z,'-','blue')
xlabel('k')
ylabel('z_k')
grid('on')
axis([0 1000 -20 20])
subplot(122)
plot(Z,Zp,'o','MarkerSize',2,'MarkerEdgeColor','blue','MarkerFaceColor',...
                             'blue')
xlabel('z_k')
ylabel('z_{k+1}')
axis([-20 20 -32 17])
grid on
hold on
pa = plot(xdom,lambdaslope*xdom,'green','LineWidth',1.3);
pb = plot(xdom,A*xdom,'red','LineWidth',1.3);
set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');