dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/";
dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtIII_1k/";

%dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWcExcl_mtI/";
%dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtII/";


addpath('/home/cs07/Scripts/FUNCTIONS_scripts');
addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT');

starting = 0;
ending   = 999;  %999;

ent_arr = zeros(ending+1, 6, 'double');

for n = starting:ending
    file_path_se        = append( dat_path_se, sprintf("SN_%d_*.mat", round(n)) );
    file_path_L1        = append( dat_path_L1, sprintf("SN_%d_*.mat", round(n)) );

    folder_info_se      = dir(file_path_se);
    folder_info_L1      = dir(file_path_L1);
    file_info_se        = append(dat_path_se, folder_info_se(1, 1).name);
    file_info_L1        = append(dat_path_L1, folder_info_L1(1, 1).name);
    Training_History_se = load(file_info_se,'Training_History').Training_History;
    Training_History_L1 = load(file_info_L1,'Training_History').Training_History;


    Aij_se         = cell2mat(Training_History_se(end, 1));
    accuracy_se    = cell2mat(Training_History_se(end, 6));
    accuracy_L1    = cell2mat(Training_History_L1(end, 6));
    Aij_L1         = cell2mat(Training_History_L1(end, 1));

    output_se  = KS_ent(Aij_se);
    output_L1  = KS_ent(Aij_L1);
    
    [q_se, R_se]   = Ricci(abs(Aij_se), 1);
    [q_L1, R_L1]   = Ricci(abs(Aij_L1), 1);

    ent_arr(n+1, 1) = output_se;
    ent_arr(n+1, 2) = q_se;
    ent_arr(n+1, 3) = output_L1;
    ent_arr(n+1, 4) = q_L1;
    ent_arr(n+1, 5) = accuracy_se;
    ent_arr(n+1, 6) = accuracy_L1;

end

ent_arr_se = ent_arr(ent_arr(:, 5) >=0.9, :); %comment out for accuracy
ent_arr_L1 = ent_arr(ent_arr(:, 6) >=0.9, :);

figure1 = scatter(ent_arr_se(:, 2), ent_arr_se(:, 1), 'filled')
hold on
title('seRNN: L = L_{task} + W \cdot D \cdot C')
%title('seRNN: L = L_{task} + W \cdot C')
xlabel('Geometric Modularity')
ylabel('Topological Entropy') %Shannon Entropy
grid on

%axis([0 1 2 3.5])
%axis([0 1 0 5])
%axis([0 0.4 0 1])

figure1.MarkerFaceColor = [243,175,182]./255;
b = gca;
b.TickDir = 'out';
b.TickLength = [.02 .02];
b.FontName = 'Arial';
b.FontSize = 17;


saveas(figure1, '/home/cs07/geoQ90_KS_se.pdf')
hold off

figure2 = scatter(ent_arr_L1(:, 4), ent_arr_L1(:, 3), 'filled')
hold on
title('L1: L = L_{task} + W')
%title('L1: L = L_{task} + W \cdot D')
xlabel('Geometric Modularity')
ylabel('Topological Entropy')
grid on

%axis([0 1 2 3.5])
%axis([0 1 0 5])
%axis([0 0.4 0 1])

figure2.MarkerFaceColor = [121,218,216]./255;
b = gca;
b.TickDir = 'out';
b.TickLength = [.02 .02];
b.FontName = 'Arial';
b.FontSize = 17;

saveas(figure2, '/home/cs07/geoQ90_KS_L1.pdf')
hold off
