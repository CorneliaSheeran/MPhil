% OR-curvature
% adj
% information matrix (syn, red)

dat_path_se = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/";
dat_path_L1 = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtIII_1k/";

addpath('/home/cs07/Scripts/FUNCTIONS_scripts');
addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT');

% starting = 0;
% ending   = 999;

networks = [10, 20, 400, 500, 890, 900];

ent_arr = zeros(ending+1, 6, 'double');

for n = networks

    file_path_se        = append( dat_path_se, sprintf("SN_%d_*.mat", round(n)) );
    file_path_L1        = append( dat_path_L1, sprintf("SN_%d_*.mat", round(n)) );

    folder_info_se      = dir(file_path_se);
    folder_info_L1      = dir(file_path_L1);
    file_info_se        = append(dat_path_se, folder_info_se(1, 1).name);
    file_info_L1        = append(dat_path_L1, folder_info_L1(1, 1).name);
    Training_History_se = load(file_info_se,'Training_History').Training_History;    
    Training_History_L1 = load(file_info_L1,'Training_History').Training_History;

    
    Aij_se       = cell2mat(Training_History_se(end, 1));
    Aij_L1       = cell2mat(Training_History_L1(end, 1));
    accuracy_se  = cell2mat(Training_History_se(end, 6));
    accuracy_L1  = cell2mat(Training_History_L1(end, 6));

    sprintf("/home/cs07/PID_eigen/se_syn_%d.csv", round(n))

    se_syn_path = sprintf("/home/cs07/PID_eigen/se_syn_%d.csv", round(n));
    se_red_path = sprintf("/home/cs07/PID_eigen/se_red_%d.csv", round(n));
    L1_syn_path = sprintf("/home/cs07/PID_eigen/L1_syn_%d.csv", round(n));
    L1_red_path = sprintf("/home/cs07/PID_eigen/L1_red_%d.csv", round(n));

    se_syn = csvread(se_syn_path);
    se_red = csvread(se_red_path);
    L1_syn = csvread(L1_syn_path);
    L1_red = csvread(L1_red_path);

    se_syn = reshape(se_syn, 100, 100, 50);
    se_red = reshape(se_red, 100, 100, 50);
    L1_syn = reshape(L1_syn, 100, 100, 50);
    L1_red = reshape(L1_red, 100, 100, 50);

    [q_se, R_se]   = Ricci(abs(Aij_se), 1);
    [q_L1, R_L1]   = Ricci(abs(Aij_L1), 1);
   
    to_eigen = [Aij_se, Aij_L1, R_se, R_L1, se_syn, se_red, L1_syn, L1_red];

    for matrixed=to_eigen
        find_eig(matrixed)
    end

    ent_arr(n+1, 1) = output_se;
    ent_arr(n+1, 2) = modu_se;
    ent_arr(n+1, 3) = output_L1;
    ent_arr(n+1, 4) = modu_L1;
    ent_arr(n+1, 5) = accuracy_se;
    ent_arr(n+1, 6) = accuracy_L1;
    
end

ent_arr_se = ent_arr(ent_arr(:, 5) >=0.9);
ent_arr_L1 = ent_arr(ent_arr(:, 6) >=0.9);


function [max_vec, largest_value] = find_eig(M)

[data_vec, data_val] = eig(M);

eig_data = zeros(N, 1);

for n=1:100
    eig_data(n, 1) = data_val(n, n);
end

[sorted_values, sorted_indices] = sort(eig_data, 'descend');
largest_value = sorted_values(1);
largest_index = sorted_indices(1);

max_vec = data_vec(:, largest_index);

end
