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

starting = 1; %1;
ending   = 99; %99

ent_arr = zeros(ending+1, 6, 'double');

for n = starting:ending

    if n <10 
        bse = n/100000; %1000=S
        cse = num2str(bse);
        bl1 = n/100000;
        cl1 = num2str(bl1);
    else
	    if mod(n,10) == 0
	        bse = n/100000; %1000=S
            cse = sprintf('%.4f', bse);
            bl1 = n/100000;
            cl1 = sprintf('%.4f', bse);
    
        else
            bse = n/100000; %1000=S
            cse = sprintf('%.5f', bse);
            bl1 = n/100000;
            cl1 = sprintf('%.5f', bse);
	    
	    end
    end

    se_acc_full  = csvread(dat_path_se_acc + 'mem_acc_' + cse + '_.csv');
    L1_acc_full  = csvread(dat_path_L1_acc + 'mem_acc_' + cl1 +'_.csv');

    file_info_se = dat_path_se + 'se_' + int2str(n) + '_epoch=50.mat';
    file_info_L1 = dat_path_se + 'L1_' + int2str(n) + '_epoch=50.mat';

    Aij_se       = load(file_info_se,'numpy_array').numpy_array;    
    Aij_L1       = load(file_info_L1,'numpy_array').numpy_array;    
    accuracy_se  = se_acc_full(end);
    accuracy_L1  = L1_acc_full(end);

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


saveas(figure1, '/home/cs07/DVS_geoQ90_KS_se.pdf')
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

saveas(figure2, '/home/cs07/DVS_geoQ90_KS_L1.pdf')
hold off
