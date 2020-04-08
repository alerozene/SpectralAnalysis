function dxdt = P4T1_fx(t,x)
dxdt = zeros(2,1);
dxdt(1) = x(1);
dxdt(2) = x(2)+x(1)^2;