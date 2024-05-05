dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/";
dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtIII_1k/";

addpath('/home/cs07/Scripts/FUNCTIONS_scripts');
addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT');

starting = 1;
ending   = 1000;
space    = 10;
tau_end  = 50;

dim_1   = round((ending-starting+1)/space)
ent_arr = zeros(dim_1, 4, round(tau_end/5), 'double');

for tau = 1:tau_end
    for n = starting:space:ending

        file_path_se        = append( dat_path_se, sprintf("SN_%d_*.mat", round(n-1)) );
        file_path_L1        = append( dat_path_L1, sprintf("SN_%d_*.mat", round(n-1)) );
    
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
        
        [q_se, R_se]   = Ricci(abs(Aij_se), tau);
        [q_L1, R_L1]   = Ricci(abs(Aij_L1), tau);
        
        pos_n = round(n/space) + 1       

        ent_arr(pos_n, 1, tau) = q_se;
        ent_arr(pos_n, 2, tau) = q_L1;
        ent_arr(pos_n, 3, tau) = accuracy_se;
        ent_arr(pos_n, 4, tau) = accuracy_L1;

    end
end

ent_arr_se = ent_arr(ent_arr(:, 3, 1) >=0.9, :, :); %comment out for accuracy
ent_arr_L1 = ent_arr(ent_arr(:, 4, 1) >=0.9, :, :);

x_tau = 1:tau_end;

ent_arr_se = mean(ent_arr_se, 1); %average over networks for accuracy > 90
size(ent_arr_se)

figure1 = scatter(x_tau, ent_arr_se(1, 1, :), 'filled')
hold on
%title('KS/Topological Entropy Against Modularity') %Shannon Entropy Against Modularity
title('seRNN: L = L_{task} + W \cdot D \cdot C')
xlabel('Tau')
ylabel('Geometric Modularity') %Shannon Entropy
grid on

figure1.MarkerFaceColor = [243,175,182]./255;
b = gca;
b.TickDir = 'out';
b.TickLength = [.02 .02];
b.FontName = 'Arial';
b.FontSize = 17;

saveas(figure1, '/home/cs07/tau90_geoQ_se.pdf')
hold off


ent_arr_L1 = average(ent_arr_L1, 1);

figure2 = scatter(x_tau, ent_arr_L1(1, 2, :), 'filled')
hold on
title('L1: L = L_{task} + W')
xlabel('Tau')
ylabel('Geometric Modularity')
grid on

figure2.MarkerFaceColor = [121,218,216]./255;
b = gca;
b.TickDir = 'out';
b.TickLength = [.02 .02];
b.FontName = 'Arial';
b.FontSize = 17;

saveas(figure2, '/home/cs07/tau90_geoQ_L1.pdf')
hold off
