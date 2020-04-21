% SC42135 Spectral Analysis of Nonlinear/Infinite-Dimensional Systems.
% Practice Session 1. Toeplitz operators
% Required: chebfun package from chebfun.org
clc; close all
%% Task 1 Matrix representation of T*Pn and compute the norm ||T*Pn||
% State space system
c = [1 -1 0]';
b = [1 0 -2]';
A = [0.5 3 2; 0 -0.5 -1; 0 0 0.2];
d = 0.2;
% Iterate to find convergence of norm
N = 100;
cnvrg = zeros(N,1);
for ii=1:N
    h=zeros(ii,1);
    h(1)=d;
    h2 = h;
    
    for ii=0:ii-1
        h(ii+2)=c'*(A^ii)*b;
    end
    
    T = toeplitz(h,h2);
    cnvrg(ii) = norm(T,2);
end
%% Task 2 Use chebfun and max to compute ||T||
H = chebfun(@(w) d+exp(-1i*w)*c'*((eye(3)-exp(1i*w)*A)^(-1))*b,...
    [0 2*pi]);
max(H, 'global')
T_normChebfun = norm(max(H, 'global'));
%% Figures
figure;
grid on
plot(1:limit,cnvrg,'r--','lineWidth',1.5)
hold on
plot(ones(limit,1)*T_normChebfun,'b-','lineWidth',1.5)
legend('Upper norm','Chebfun')
