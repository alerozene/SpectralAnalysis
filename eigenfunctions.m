function [gammadelta,psivec] = eigenfunctions(a,u,lambda)
% Given a signal u and eigenvalue(s), recover eigenfunction and exponential 
% coefficient at point a 

for ii = 1:length(lambda)
    
    N = chebop(-a,a);
    % Schrödinger operator
    N.op = @(z,psi) diff(psi,2) - (lambda(ii)-u)*psi;
    % Maybe this was a struggle (changing the i.c: i though -a = 1 was sufficiently large)
    N.lbc = [exp(-1i*sqrt(-lambda(ii))*-a); ...
             -1i*sqrt(-lambda(ii))*exp(-1i*sqrt(-lambda(ii))*-a)];
    % solve task 1 chebfun 10.2 style
    psi = N\0;
    delta_psi = diff(psi);
    psivec(:,ii) = [psi(a);delta_psi(a)];
    
    % eigenfunction as linear combination of exponentials (task 2)
    A = zeros(2);
    A(1,1)= exp(1i*sqrt(-lambda(ii))*a);
    A(1,2)= exp(-1i*sqrt(-lambda(ii))*a);
    A(2,1)= 1i*sqrt(-lambda(ii))*exp(1i*sqrt(-lambda(ii))*a);
    A(2,2)= -1i*sqrt(-lambda(ii))*exp(-1i*sqrt(-lambda(ii))*a);
    gammavector = A\[psi(a); delta_psi(a)];
    gammadelta(:,ii) = gammavector;  
    
end
