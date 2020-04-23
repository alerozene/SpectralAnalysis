% SC42135 Spectral Analysis of Nonlinear/Infinite-Dimensional Systems.
% Final assignment
close all; clear all; clc
%% Task 1. IVP to recover eigenfunction using chebfun package[Section 10.2]

a = 2;
lambda = -1;
u = 1;

N = chebop(-a,a);
N.op = @(psi) diff(psi,2) - (lambda-u)*psi;
init_cndtn = [exp(-1i*sqrt(-lambda)*-a); -1i*sqrt(-lambda)*exp(-1i*sqrt(-lambda)*-a)];
N.lbc = init_cndtn;

psi = N\0;
delta_psi = diff(psi);
psi(a)
delta_psi(a)

%% Task 2. 

A = zeros(2);
A(1,1)= exp(1i*sqrt(-lambda)*a);
A(1,2)= exp(-1i*sqrt(-lambda)*a);
A(2,1)= 1i*sqrt(-lambda)*exp(1i*sqrt(-lambda)*a);
A(2,2)= -1i*sqrt(-lambda)*exp(-1i*sqrt(-lambda)*a);

gammadelta = A\[psi(a); delta_psi(a)];

%% Task 3.
tic
aa = 1;
a = 25;
for lambda = -25:0.1:0
    N = chebop(-a,a);
    d = [-a,a];
    z = chebfun('z', d);
    u = 20*exp(-z^4);
    N.op = @(z,psi) diff(psi,2) - (lambda-u)*psi;
    init_cndtn = [exp(-1i*sqrt(-lambda)*-a); -1i*sqrt(-lambda)*exp(-1i*sqrt(-lambda)*-a)];
    N.lbc = init_cndtn;
    psi = N\0;
    delta_psi = diff(psi);
    A = zeros(2);
    A(1,1)= exp(1i*sqrt(-lambda)*a);
    A(1,2)= exp(-1i*sqrt(-lambda)*a);
    A(2,1)= 1i*sqrt(-lambda)*exp(1i*sqrt(-lambda)*a);
    A(2,2)= -1i*sqrt(-lambda)*exp(-1i*sqrt(-lambda)*a);
    gammadelta = A\[psi(a); delta_psi(a)];
    
    gd(:,aa) = gammadelta;
    aa = aa+1;
    
end
toc

%% Task 4.Write a Matlab script that recreates Figure 1 
% data from 3.2.1 p11
lambda = 0.5;
n_variance = sqrt(10);
k = 1e3;
nk = n_variance*randn([1 k]);
z = zeros([1 k]);
% eq(18)
for ii = 2:k
    z(ii) = lambda*z(ii-1) +nk(ii);
end

% Domains for Z and Zp
dom1 = 1:k-1;
dom2 = 2:k;
% Data amtrices
Z = z(1:k-1);
Zp = z(2:k);


% 2.2
X = Z;
Y = Zp;
% Definition 1 (lamda hat) 
A = Y*pinv(X);

figure(4)
subplot(121)
plot(z,'-')
xlabel('k')
ylabel('z_k')
grid('on')
axis([0 1000 -20 20])
subplot(122)
plot(Z,Zp,'o')
xlabel('z_k')
ylabel('z_{k+1}')
axis([-20 20 -32 17])
hold on
plot(zhat,'-')


