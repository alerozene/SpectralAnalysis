function Y = P4T4(X,d)
% Compute measurements x1^p*x2^q for a matrix of states


% Allocate according to Task 4 in practice session sheet
for kk = 1:size(X,2)
    a = 1;
    for ii =0:d
        for jj=0:d
            if ii+jj<=d
                Y(a,kk)= X(1,kk)^(ii)*X(2,kk)^(jj);
                a = a+1;
            end
        end
    end
end
