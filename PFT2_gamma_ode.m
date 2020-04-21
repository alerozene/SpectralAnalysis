function dxda = PFT2_gamma_ode(a,psi)
A(1,1)= exp(1i*sqrt(-lambda)*a);
A(1,2)= exp(-1i*sqrt(-lambda)*a);
A(2,1)= 1i*sqrt(-lambda)*exp(1i*sqrt(-lambda)*a);
A(2,2)= -1i*sqrt(-lambda)*exp(-1i*sqrt(-lambda)*a);
dxda = A*psi;
