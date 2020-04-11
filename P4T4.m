function Y = P4T4(X,d)
% Compute measurements x1^p*x2^q for a matrix of states


% Allocate according to Task 4 in practice session sheet
for kk = 1:length(X)
    a = 1;
    for ii =0:d
        for jj=0:d
            if ii+jj<=d
                Y(kk,a)= X(kk,1)^(ii)*X(kk,2)^(jj);
                a = a+1;
            end
        end
    end
end
