function communicability = CalCom(A)
    % Calculate weighted network communicability
    A  = abs(A);
    D  = zeros(N);
    dd = sum(A, 2);

    for m = 1:N
        D(m, m) = dd(m);
    end

    M = sqrtm(inv(D))*A*sqrtm(inv(D));

    communicability = expm(M);
end