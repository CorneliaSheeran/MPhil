%Modularity Against Entropy for split loss-function

dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWcExcl_mtI/";
dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtII/";

%mazeGenI_SE1_sWcExcl_mtI
%mazeGenI_L1_mtII

addpath('/home/cs07/Scripts/FUNCTIONS_scripts');
addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT');

starting = 0;
ending   = 99;

ent_arr = zeros(ending+1, 6, 'double');

for n = starting:ending
    
   %file_path        = append( dat_path, sprintf("SN_%d_*.mat", round(n)) );
    %folder_info      = dir(file_path);
    %file_info        = append(dat_path, folder_info(1, 1).name);
    %Training_History = load(file_info,'Training_History').Training_History; 
    %Reg_strength     = load(file_info,'Regulariser_Strength').Regulariser_Strength; 
    


        file_path_se        = append( dat_path_se, sprintf("SN_%d_*.mat", round(n)) );
    file_path_L1        = append( dat_path_L1, sprintf("SN_%d_*.mat", round(n)) );

folder_info_se      = dir(file_path_se);
folder_info_L1      = dir(file_path_L1);
file_info_se        = append(dat_path_se, folder_info_se(1, 1).name);
file_info_L1        = append(dat_path_L1, folder_info_L1(1, 1).name);
    Training_History_se = load(file_info_se,'Training_History').Training_History;    
    Training_History_L1 = load(file_info_L1,'Training_History').Training_History;

    
    Aij_se     = cell2mat(Training_History_se(end, 1));
accuracy_se    = cell2mat(Training_History_se(end, 6));
accuracy_L1    = cell2mat(Training_History_L1(end, 6));
    Aij_L1     = cell2mat(Training_History_L1(end, 1));

output_se  = RW_ent(Aij_se);
output_L1  = RW_ent(Aij_L1);

    [C_se, modu_se]    = modularity_und(abs(Aij_se));
[C_L1, modu_L1]    = modularity_und(abs(Aij_L1));    

    ent_arr(n+1, 1) = output_se;
    ent_arr(n+1, 2) = modu_se;
    ent_arr(n+1, 3) = output_L1;
    ent_arr(n+1, 4) = modu_L1;
    ent_arr(n+1, 5) = accuracy_se;
    ent_arr(n+1, 6) = accuracy_L1;
    
end

%ent_arr = ent_arr(ent_arr(:,2) >= -1, :);
%ent_arr = ent_arr(ent_arr(:, 4) >= -1, :);

ent_arr_se = ent_arr(ent_arr(:, 5) >=0.9, :);
ent_arr_L1 = ent_arr(ent_arr(:, 6) >=0.9, :);

figure1 = scatter(ent_arr_se(:, 2), ent_arr_se(:, 1), 'filled') % [], ent_arr_se(:, 5), 'filled');
hold on
%colorbar
%scatter(ent_arr(:, 4), ent_arr(:, 3), [], ent_arr(:, 5));
title('C_{ij} only: L = L_{task} + W \cdot C') %KS/Topological Entropy Against Modularity 
xlabel('Modularity')
ylabel('Shannon Entropy')
%legend("Spatially Embedded", "L1")
grid on
axis([0 1 0 5])

figure1.MarkerFaceColor = [0 0 0]; %[243,175,182]./255;

b = gca; 
b.TickDir = 'out'; 
b.TickLength = [.02 .02]; 
b.FontName = 'Arial'; 
b.FontSize = 20;

saveas(figure1, '/home/cs07/mod_90_com_RW_plot.pdf')
hold off

figure2 = scatter(ent_arr_L1(:, 4), ent_arr_L1(:, 3), 'filled') % [], ent_arr_L1(:, 6), 'filled');
hold on
%colorbar
%scatter(ent_arr(:, 4), ent_arr(:, 3), [], ent_arr(:, 5));
title('D_{ij} only: L = L_{task} + W \cdot D') %KS/Topological Entropy Against Modularity 
xlabel('Modularity')
ylabel('Shannon Entropy')
%legend("Spatially Embedded", "L1")
grid on

axis([0 1 0 5])

figure2.MarkerFaceColor = [0.4660 0.6740 0.1880];  %[121,218,216]./255;

b = gca;
b.TickDir = 'out';
b.TickLength = [.02 .02]; 
b.FontName = 'Arial';
b.FontSize = 20;

saveas(figure2, '/home/cs07/mod_90_spatial_RW_plot.pdf')
hold off
