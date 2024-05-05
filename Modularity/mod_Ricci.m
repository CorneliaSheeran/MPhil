dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/";
dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtIII_1k/";

%dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWcExcl_mtI/";
%dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtII/";


addpath('/home/cs07/Scripts/FUNCTIONS_scripts');
addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT');

starting = 0;
ending   = 999;  %999;

ent_arr = zeros(ending+1, 6, 100, 'double');
N = 100;
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

    %output_se  = RW_ent(Aij_se);
    %output_L1  = RW_ent(Aij_L1);

    [C_se, modu_se]    = modularity_und(abs(Aij_se));
    [C_L1, modu_L1]    = modularity_und(abs(Aij_L1));
    

    [q_se, R_se]   = Ricci(Aij_se, 1);
    [q_L1, R_L1]   = Ricci(Aij_L1, 1);
    R_se = sum(abs(R_se), 1)/N;
    R_L1 = sum(abs(R_L1), 1)/N;

    ent_arr(n+1, 1, :) = modu_se*ones(1, N);
    ent_arr(n+1, 2, :) = R_se;
    ent_arr(n+1, 3, :) = modu_L1*ones(1, N);
    ent_arr(n+1, 4, :) = R_L1;
    ent_arr(n+1, 5, :) = accuracy_se*ones(1, N);
    ent_arr(n+1, 6, :) = accuracy_L1*ones(1, N);

end

ent_arr_se = ent_arr(ent_arr(:, 5, 1) >=0.9, :, :); %comment out for accuracy
ent_arr_L1 = ent_arr(ent_arr(:, 6, 1) >=0.9, :, :);

ent_arr_se = ent_arr_se(:, 1:2, :);
[~,index_se] = sort(ent_arr_se(:, 1, 1));
ent_arr_se = ent_arr_se(index_se, 2, :); 

ent_arr_L1 = ent_arr_L1(:, 1:2, :);
[~,index_L1] = sort(ent_arr_L1(:, 1, 1));
ent_arr_L1 = ent_arr_L1(index_L1, 2, :);

indexed = 1:10:1000;
colour_data = reshape(ent_arr_se, 100, []);
colour_data = colour_data(:, indexed);

figure1 = heatmap(colour_data, 'Colormap', parula)
xlabel('Increasing Modularity');
ylabel('Nodes');
b = gca;
b.FontName = 'Arial';
b.FontSize = 17;
saveas(figure1, '/home/cs07/mod90_absR3D_se.pdf')

colour_data2 = reshape(ent_arr_L1, 100, []);
colour_data2 = colour_data2(:, indexed);

figure2 = heatmap(colour_data2, 'Colormap', parula)
xlabel('Increasing Modularity');
ylabel('Nodes');
b = gca;
b.FontName = 'Arial';
b.FontSize = 17;
saveas(figure2, '/home/cs07/mod90_absR3D_L1.pdf')
















for dont_go_here = []

[X,Y] = meshgrid(-1:1, 1:100);
Z = ent_arr_se(:, 2, :);
figure1 = surf(X,Y,Z)
hold on
colorbar
saveas(figure1, '/home/cs07/Riccimod_3D_se.pdf')
hold off

X = 1:100;
Y = ent_arr_L1(:, 2, :)
Z = ent_arr_L1(:, 1, :);
C = ent_arr_L1(:, 1, :);
figure1 = surf(X,Y,Z,C)
hold on
colorbar
saveas(figure1, '/home/cs07/Riccimod_3D_L1.pdf')
hold off

figure1 = plot(ent_arr_se(:, 2, 1), ent_arr_se(:, 1, 1))
hold on
for i=2:10
plot(ent_arr_se(:, 2, i), ent_arr_se(:, 1, i))
end
title('seRNN: L = L_{task} + W \cdot D \cdot C')
%title('seRNN: L = L_{task} + W \cdot C')
xlabel('Ricci Curvature')
ylabel('Modularity') %Shannon Entropy
grid on

%axis([0 1 2 3.5])
%axis([0 1 0 5])
%axis([0 0.4 0 1])

figure1.MarkerFaceColor = [243,175,182]./255;
%b = gca;
%b.TickDir = 'out';
%b.TickLength = [.02 .02];
%b.FontName = 'Arial';
%b.FontSize = 17;


saveas(figure1, '/home/cs07/Riccimod_90_se.pdf')
hold off

figure2 = plot(ent_arr_L1(:, 4, 1), ent_arr_L1(:, 3, 1))
hold on
for i=2:10
plot(ent_arr_se(:, 4, i), ent_arr_se(:, 3, i))
end
%title('KS/Topological Entropy Against Modularity') %Shannon Entropy Against Modularity
title('L1: L = L_{task} + W')
%title('L1: L = L_{task} + W \cdot D')
xlabel('Ricci Curvature')
ylabel('Modularity')
grid on

%axis([0 1 2 3.5])
%axis([0 1 0 5])
%axis([0 0.4 0 1])

figure2.MarkerFaceColor = [121,218,216]./255;
%b = gca;
%b.TickDir = 'out';
%b.TickLength = [.02 .02];
%b.FontName = 'Arial';
%b.FontSize = 17;

saveas(figure2, '/home/cs07/Riccimod_90_L1.pdf')
hold off
end 
