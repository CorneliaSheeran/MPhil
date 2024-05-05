%compute the Tononi, Sporns and Edeleman (TSE) Entropy

function H = tse_com(X, min_k, max_k)
    %input: X = matrix of time series 
    %ouput: H = TSE entropy

    n      = size(X, 1);
    mask   = 1:1:n;
    CN_all = zeros((max_k-min_k)+1, 1);

    Y      = reshape(X, size(X, 1)*size(X, 2), 1);
    size(Y)
    Y(1:5, 1)   
 H_tot  = shan_ent(Y);

    for k=min_k:max_k
        %partion of size k:
        index = nchoosek(mask, k);
        iter  = size(index, 1);
        av    = zeros(iter, 1);
    
        for i=1:iter
            part   = index(i, :);
            dat_bi = X(part, :);
            
            ent   = shan_ent(reshape(dat_bi, 1, size(dat_bi, 1)*size(dat_bi, 2)));
            av(i) = ent;

        end

        H_av   = mean(av);
        CN_all(k, 1) = H_av - (k/n)*H_tot;
        
    end

    H = sum(CN_all);
end
