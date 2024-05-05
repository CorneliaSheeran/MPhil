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
epochs_n = 50; %50

ent_arr = zeros(ending+1, epochs_n, 4, 'double'); %50 epochs

for n = starting:ending

    if n <10 
        bse = n/100000; %1000=S
        cse = num2str(bse);
        bl1 = n/100000
        cl1 = num2str(bl1)
    else
	if mod(n,10) == 0
	
	bse = n/100000; %1000=S
        cse = sprintf('%.4f', bse);
        bl1 = n/100000
        cl1 = sprintf('%.4f', bse)

	else

        bse = n/100000; %1000=S
        cse = sprintf('%.5f', bse);
        bl1 = n/100000
        cl1 = sprintf('%.5f', bse)
	
	end
    end

    se_acc_full = csvread(dat_path_se_acc + 'mem_acc_' + cse + '_.csv');
    L1_acc_full = csvread(dat_path_L1_acc + 'mem_acc_' + cl1 +'_.csv');
    
    for ep = 1:epochs_n
        accuracy_se = se_acc_full(ep);
        accuracy_L1 = L1_acc_full(ep);       

        file_info_se        = dat_path_se + 'se_' + int2str(n) + '_epoch=' + int2str(ep);
        file_info_L1        = dat_path_se + 'L1_' + int2str(n) + '_epoch=' + int2str(ep);
        Training_History_se = load(file_info_se,'numpy_array').numpy_array;    
        Training_History_L1 = load(file_info_L1,'numpy_array').numpy_array;

        Aij_se      = Training_History_se;
        Aij_L1      = Training_History_L1;
                
       % [q_se, R_se]   = Ricci(Aij_se, 1);
       % [q_L1, R_L1]   = Ricci(Aij_L1, 1);

        output_se  = RW_ent(Aij_se);
        output_L1  = RW_ent(Aij_L1);
    
        ent_arr(n+1, ep, 1) = output_se;
        ent_arr(n+1, ep, 2) = output_L1;
        ent_arr(n+1, ep, 3) = accuracy_se;
        ent_arr(n+1, ep, 4) = accuracy_L1;
    end 
    
end

ent_arr_se = ent_arr(ent_arr(:, end, 3) >=0.9, :, :);
ent_arr_L1 = ent_arr(ent_arr(:, end, 3) >=0.9, :, :);

epochs = 1:epochs_n;

ent_arr_se = ent_arr_se(:, :, 1);
ent_arr_L1 = ent_arr_L1(:, :, 2);

error_se   = reshape(std(ent_arr_se), epochs_n, 1);
error_L1   = reshape(std(ent_arr_L1), epochs_n, 1);
ent_arr_se = reshape(mean(ent_arr_se), epochs_n, 1);
ent_arr_L1 = reshape(mean(ent_arr_L1), epochs_n, 1);

figure1 = errorbar(epochs, ent_arr_se, error_se, 'Color', [243,175,182]./255)
hold on
title('seRNN: L = L_{task} + W \cdot D \cdot C')
xlabel('Epochs')
ylabel('Shannon Entropy')
grid on

%figure1.MarkerFaceColor = [243,175,182]./255;
b = gca;
b.TickDir = 'out'; 
b.TickLength = [.02 .02]; 
b.FontName = 'Arial'; 
b.FontSize = 17;

saveas(figure1, '/home/cs07/DVS_epoch90_RW_se.pdf')
hold off


figure2 = errorbar(epochs, ent_arr_L1, error_L1, 'Color', [121,218,216]./255)
hold on
title('L1: L = L_{task} + W')  
xlabel('Epochs')
ylabel('Shannon Entropy')
grid on

%figure2.MarkerFaceColor = [121,218,216]./255;
b = gca;
b.TickDir = 'out'; 
b.TickLength = [.02 .02]; 
b.FontName = 'Arial'; 
b.FontSize = 17;

saveas(figure2, '/home/cs07/DVS_epoch90_RW_L1.pdf')
hold off
