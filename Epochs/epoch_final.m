dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/";
dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtIII_1k/";

%mazeGenI_SE1_sWcExcl_mtI
%mazeGenI_L1_mtII

addpath('/home/cs07/Scripts/FUNCTIONS_scripts');
addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT');

starting = 0;
ending   = 999;

ent_arr = zeros(ending+1, 4, 100, 10, 'double'); %10 epochs

for n = starting:ending

    file_path_se        = append( dat_path_se, sprintf("SN_%d_*.mat", round(n)) );
    file_path_L1        = append( dat_path_L1, sprintf("SN_%d_*.mat", round(n)) );
    folder_info_se      = dir(file_path_se);
    folder_info_L1      = dir(file_path_L1);
    file_info_se        = append(dat_path_se, folder_info_se(1, 1).name);
    file_info_L1        = append(dat_path_L1, folder_info_L1(1, 1).name);
    Training_History_se = load(file_info_se,'Training_History').Training_History;    
    Training_History_L1 = load(file_info_L1,'Training_History').Training_History;

    for e = 1:10
        Aij_se      = cell2mat(Training_History_se(e, 1));
        accuracy_se = cell2mat(Training_History_se(e, 6));
        accuracy_L1 = cell2mat(Training_History_L1(e, 6));
        Aij_L1      = cell2mat(Training_History_L1(e, 1));
    
        output_se  = RW_ent(Aij_se);
        output_L1  = RW_ent(Aij_L1);
    
        ent_arr(n+1, 1, e) = output_se;
        ent_arr(n+1, 2, e) = output_L1;
        ent_arr(n+1, 3, e) = accuracy_se;
        ent_arr(n+1, 4, e) = accuracy_L1;
    end 
    
end

%ent_arr_se = ent_arr(ent_arr(:, 5) >=0.9, :);
%ent_arr_L1 = ent_arr(ent_arr(:, 6) >=0.9, :);

epochs = 1:10;
ent_arr_se = average(ent_arr(:, 1, :), 1);
ent_arr_L1 = average(ent_arr(:, 2, :), 1);

error_se = std(ent_arr(:, 1, :), 0, 1);
error_L1 = std(ent_arr(:, 2, :), 0, 1);

figure1 = errorbar(epochs, ent_arr_se, error_se)
hold on
title('seRNN: L = L_{task} + W \cdot D \cdot C')
xlabel('Epochs')
ylabel('Shannon Entropy')
grid on

figure1.MarkerFaceColor = [243,175,182]./255;
%b = gca;
%b.TickDir = 'out'; 
%b.TickLength = [.02 .02]; 
%b.FontName = 'Arial'; 
%b.FontSize = 17;

saveas(figure1, '/home/cs07/mod_RW_se.pdf')
hold off


figure1 = errorbar(epochs, ent_arr_L1, error_L1)
hold on
title('L1: L = L_{task} + W')  
xlabel('Epochs')
ylabel('Shannon Entropy')
grid on

figure2.MarkerFaceColor = [121,218,216]./255;%b = gca;
%b.TickDir = 'out'; 
%b.TickLength = [.02 .02]; 
%b.FontName = 'Arial'; 
%b.FontSize = 17;

saveas(figure2, '/home/cs07/mod_RW_L1.pdf')
hold off