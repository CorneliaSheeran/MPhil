function [expectedC,expectedL] = ER_Expected_L_C(k,n,nperm,m)
%   
% ER_EXPECTED_L_C the expected path-length and clustering of an ER random graph
%
% [C,L] = ER_EXPECTED_L_C(K,n,nperm) for a network of n nodes and mean degree K, 
% computes the expected shortest path-length L and expected clustering 
% coefficient C if that network was an Erdos-Renyi random graph
% 
% K = mean degree
%
% n = number of nodes in the network
%
% m = number of the edges in the network
% 
% Mark Humphries 3/2/2017

% initalise
clu_perm = zeros(1,nperm);
cpl_perm = zeros(1,nperm);
for perm = 1:nperm;
    % make the ER graph
    Arand = makerandCIJ_und(n,m);
    % compute the CC
    clu_perm(perm) = mean(clustering_coef_bu(Arand));
    % compute the CPL
    cpl_perm(perm) = charpath(Arand);
end
expectedC = mean(clu_perm);
expectedL = mean(cpl_perm);

%{
expectedC = k/n;
z1 = k;
z2 = k^2;
expectedL = (log((n-1)*(z2 - z1) + z1^2) - log(z1^2)) / log(z2/z1);
%}