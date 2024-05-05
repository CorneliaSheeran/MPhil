%dat_path = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/";
dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/";
dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtIII_1k/";

addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT') 
addpath('/home/cs07/Scripts/FUNCTIONS_scripts')

starting = 0;
ending   = 999;

ent_arr = zeros(ending+1, 6, 'double');

for n = starting:ending


    file_path_se        = append( dat_path_se, sprintf("SN_%d_*.mat", round(n)) );
    file_path_L1        = append( dat_path_L1, sprintf("SN_%d_*.mat", round(n)) );

    folder_info_se      = dir(file_path_se);
    folder_info_L1      = dir(file_path_L1);
    file_info_se        = append(dat_path_se, folder_info_se(1, 1).name);
    file_info_L1        = append(dat_path_L1, folder_info_L1(1, 1).name);
    
    Training_History_se = load(file_info_se,'Training_History').Training_History;
    Reg_strength_se     = load(file_info_se,'Regulariser_Strength').Regulariser_Strength;
    Training_History_L1 = load(file_info_L1,'Training_History').Training_History;
    Reg_strength_L1     = load(file_info_L1,'Regulariser_Strength').Regulariser_Strength;
    
  
    Aij_se     = cell2mat(Training_History_se(end, 1));
    Aij_L1     = cell2mat(Training_History_L1(end, 1));
    
    %[q_se, R_se]   = Ricci(Aij_se, 1);
    %[q_L1, R_L1]   = Ricci(Aij_L1, 1);

    output_se  = KS_ent(Aij_se);
    output_L1  = KS_ent(Aij_L1);

    accuracy_se    = cell2mat(Training_History_se(end, 6));
    accuracy_L1    = cell2mat(Training_History_L1(end, 6));
   
    % set number of permutations for the smw null model
    nperm = 1000;
    % set the number of edges in the empirical
    % m = nnz(Aij_abs_thr)/2;
    nn = size(Aij_se,1);

    bin_A_se = threshold_proportional(abs(Aij_se), 0.1);   
    bin_A_se(bin_A_se > 0) = 1;
    m_se = nnz(bin_A_se)/2;
   
    bin_A_L1 = threshold_proportional(abs(Aij_L1), 0.1);
    bin_A_L1(bin_A_L1 > 0) = 1;
    m_L1 = nnz(bin_A_L1)/2;

    k_se = sum(bin_A_se);
    K_se = mean(k_se); % mean degree of network

    k_L1 = sum(bin_A_L1);
    K_L1 = mean(k_L1); % mean degree of network    


    [expectedC_se,expectedL_se] = ER_Expected_L_C(K_se,nn, nperm, m_se);  % L_rand and C_rand
    [expectedC_L1,expectedL_L1] = ER_Expected_L_C(K_L1,nn, nperm, m_L1);  % L_rand and C_rand
   
    [S_ws_se,C_ws_se,L_se] = small_world_ness(bin_A_se, expectedL_se, expectedC_se, 1);  % Using WS clustering coefficient
    [S_ws_L1,C_ws_L1,L_L1] = small_world_ness(bin_A_L1, expectedL_L1, expectedC_L1, 1);  % Using WS clustering coefficient
    
    ent_arr(n+1, 1) = output_se;
    ent_arr(n+1, 2) = S_ws_se;
    ent_arr(n+1, 3) = output_L1;
    ent_arr(n+1, 4) = S_ws_L1;
    ent_arr(n+1, 5) = accuracy_se;
    ent_arr(n+1, 6) = accuracy_L1;    
end

ent_arr_se = ent_arr(ent_arr(:, 5) >=0.9, :);
ent_arr_L1 = ent_arr(ent_arr(:, 6) >=0.9, :);

figure1 = scatter(ent_arr_se(:, 2), ent_arr_se(:, 1), 'filled')
hold on
title('seRNN: L = L_{task} + W \cdot D \cdot C') %title('Geometric Modularity Against Small-Worldness')
xlabel('Small-Worldness')
ylabel('Topological Entropy')
b = gca;
b.TickDir = 'out';
b.TickLength = [.02 .02];
b.FontName = 'Arial';
b.FontSize = 17;
grid on

saveas(figure1, '/home/cs07/sw90_KS_se.pdf')

hold off

figure2 = scatter(ent_arr_L1(:, 4), ent_arr_L1(:, 3), 'filled') % [], ent_arr_L1(:, 6), 'filled');
hold on
%colorbar
%scatter(ent_arr(:, 4), ent_arr(:, 3), [], ent_arr(:, 6));
title('L1: L = L_{task} + W') %title('Geometric Modularity Against Small-worldness')
xlabel('Small-Worldness')
ylabel('Topological Entropy')
b = gca;
b.TickDir = 'out';
b.TickLength = [.02 .02];
b.FontName = 'Arial';
b.FontSize = 17;
grid on

saveas(figure2, '/home/cs07/sw90_KS_L1.pdf')

hold off
