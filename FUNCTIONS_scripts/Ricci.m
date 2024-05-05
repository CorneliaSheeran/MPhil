function [Q, R] = Ricci(A, tau)
N = 100;
D = zeros(N);

%for i = 1:N
%    count = 0;
%for j = 1:N
%    if A(i,j) > 0
%       count = count + 1;
%    end
%end
%D(i,i) = count;
%end

dd = sum(A, 2);

for m = 1:N
    D(m, m) = dd(m);
end
        
L = D - A;

M = sqrtm(inv(D))*L*sqrtm(inv(D)); %This is the communicability if central object=A 

interim = expm(-M*tau);

p_list = {};

R = zeros(N);

for i = 1:N
b = interim(:,i);
p_list{i} = b;
end

for k = 1:N
for l = 1:N
    if A(k,l) > 0
        R(k,l) = 1 - (ws_distance(p_list{k}, p_list{l}, 1)/A(k,l));
    end
end
end

R(isnan(R)) = 0;

for i = 1:N
for j = 1:N
    if R(i,j) < -100
        R(i,j) = 0;
    end
end
end

maxR = max(max(R));

norm = sum(abs(R(:)));
Q = sum(R(:) - maxR)/(2*norm);
%Q = sum(R(:))/(2*norm);

end
