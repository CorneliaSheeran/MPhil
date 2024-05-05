%(DVS Gesture task):
dat_path_se = "/home/cs07/DVS/";
dat_path_L1 = "/home/cs07/DVS/";
dat_path_se_acc = "/imaging/shared/users/ja02/CBUActors/MPhil/Andrew/wrapper_seRSNN/trained_networks/seRSNN/seRSNN_DVS_run_2023_04_26/mem_acc/";
dat_path_L1_acc = "/imaging/shared/users/ja02/CBUActors/MPhil/Andrew/wrapper_seRSNN/trained_networks/L1/L1_SHD_run_2023_04_25/mem_acc/";

%(SHD task): (se, L1)
%dat_path_se = "/home/cs07/SHD/";
%dat_path_L1 = "/home/cs07/SHD/";
%dat_path_se_acc = "/imaging/shared/users/ja02/CBUActors/MPhil/Andrew/wrapper_seRSNN/trained_networks/seRSNN/seRSNN_SHD_run_2023_04_24/mem_acc/";
%dat_path_L1_acc = "/imaging/shared/users/ja02/CBUActors/MPhil/Andrew/wrapper_seRSNN/trained_networks/L1/L1_SHD_run_2023_04_25/mem_acc/";

addpath('/home/cs07/Scripts/FUNCTIONS_scripts');
addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT');

networks = [1, 25, 75, 100];

for n = networks
    file_info_se = dat_path_se + 'se_' + int2str(n)  + '_epoch=50.mat';
    file_info_L1 = dat_path_se + 'L1_' + int2str(n)  + '_epoch=50.mat';

    Aij_se	 = load(file_info_se,'numpy_array').numpy_array;
    Aij_L1	 = load(file_info_L1,'numpy_array').numpy_array;
    
    [q_se, R_se]   = Ricci(Aij_se, 1);
    [q_L1, R_L1]   = Ricci(Aij_L1, 1);

    figure1 = heatmap(R_se, 'Colormap', parula)
    xlabel('Node');
    ylabel('Node');
    b = gca;
    b.FontName = 'Arial';
    b.FontSize = 17;
    saveas(figure1, ['/home/cs07/DVS_se_Ricciheat', num2str(n), '_.pdf']);

    figure2 = heatmap(R_L1, 'Colormap', parula)
    xlabel('Node');
    ylabel('Node');
    b = gca;
    b.FontName = 'Arial';
    b.FontSize = 17;
    saveas(figure2, ['/home/cs07/DVS_L1_Ricciheat', num2str(n), '_.pdf']);

end


