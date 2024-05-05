%Compute the Shannon/Random-Walker edge entropy 
function H = SE_mat(A)
   
    addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT')
    N    = size(A, 1);
    prob = zeros(N, N);
    mu   = zeros(N, 1);
    A = abs(A);

    for i=1:N 
        for j=1:N 
            %prob(i, j) = exp(-abs(A(i, j)))/(sum(exp(-abs(A(i, :)))));
	prob(i, j) = abs(A(i, j))/(sum(abs(A(i, :))));
   
 end 
    end
    
    for m=1:N
        mu(m, 1) = - sum( prob(m, :) .* log(prob(m, :)) );
    end 
    
    H = mu;
end 
