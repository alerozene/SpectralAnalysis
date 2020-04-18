function dpsidt = P3T4(t,psi,eigvals,innervec)

u_0_reco = 0;
for ii=1:length(eigvals)
    u_0_reco = u_0_reco + eigvals(ii)*innervec(ii)*psi(2*ii);
end

for ii=1:length(eigvals)
    A(2*ii  , 2*ii-1) = 1;
    A(2*ii-1, 2*ii) = eigvals(ii)-u_0_reco;
end
dpsidt = A*psi;
