function Y = P4T4_fx(X,d)
% X is a state matrix s.t. X(dimensions,measurements at n*tau)
% d is a sclar s.t p+q<d

dxdt = [x(1); x(2)+x(1)^2];