%compute the Random Walker Kolmogorov–Sinai Network (Topology) Entropy

function H = KS_ent(A)
    addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT')
    %input: A = adjacency matrix
    %ouput: H = network entropy
    
    B      = threshold_proportional(A, 0.1);
    B(B > 0) = 1;
    lambda     = eig(B);
    lambda_mod = abs(lambda);
    [~, idx]   = max(lambda_mod);
    lambda_max = lambda(idx);
    
    H     = log(lambda_max);
end

    
