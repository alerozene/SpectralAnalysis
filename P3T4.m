function dydt = P3T4(z,y,eigvals,innervec)

for ii=1:length(eigvals)
    A(2*ii,2*ii-1) = 1;
    A(2*ii-1,2*ii) = eigvals(ii)-innervec(ii)*y(ii*2);
end
dydt = A*y;
end
