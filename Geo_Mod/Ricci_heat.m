dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/";
dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtIII_1k/";

%dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWcExcl_mtI/";
%dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtII/";


addpath('/home/cs07/Scripts/FUNCTIONS_scripts');
addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT');

networks = [0, 250, 750, 999];

%starting = 0;
%ending   = 999;  %999;

for n = networks
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

    [q_se, R_se]   = Ricci(Aij_se, 1);
    [q_L1, R_L1]   = Ricci(Aij_L1, 1);


    figure1 = heatmap(R_se, 'Colormap', parula)
    xlabel('Node');
    ylabel('Node');
    b = gca;
    b.FontName = 'Arial';
    b.FontSize = 17;
    saveas(figure1, ['/home/cs07/se_Ricciheat', num2str(n), '_.pdf']);

    figure2 = heatmap(R_L1, 'Colormap', parula)
    xlabel('Node');
    ylabel('Node');
    b = gca;
    b.FontName = 'Arial';
    b.FontSize = 17;
    saveas(figure2, ['/home/cs07/L1_Ricciheat', num2str(n), '_.pdf']);

end


