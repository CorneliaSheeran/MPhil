%function epoch_com(task_id)

dat_path = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/Activity_Recordings/";

addpath("/home/cs07/Scripts/FUNCTIONS_scripts")
addpath('/home/cs07/Scripts/BCT/2019_03_03_BCT');

multi = 1; %100

starting = 0; %round((task_id - 1)*multi);
ending   = 1; %round(task_id*multi +  1);

for n = starting:(ending)
    
    file_path   = append( dat_path, sprintf("SN_%d_*.mat", round(n)) );
    folder_info = dir(file_path);
    file_info   = append(dat_path, folder_info(1, 1).name);
    rec_data    = load(file_info,'recorded_data').recorded_data;     
    post_data_1 = rec_data(:, 8);
    post_data   = post_data_1(2:end, 1);

    num_trials  = cell2mat(rec_data(end, 1));
    tse_all     = zeros(num_trials+1, 1);

    M     = size(post_data, 1);
    width = round(M/(num_trials+1));
    length = size(cell2mat(post_data(5, 1)), 2);     

    for i=0:num_trials
        
        data_store = zeros(width, length);
        
        k = 0;
         
        for j= ((i*width)+1):(((i+1)*width))
            k = k + 1;
            data_store(k, :) = cell2mat(post_data(j, 1));
        end
                 
        tse_all(i+1, 1) = tse_com(data_store, 1, width);
    end
    
    figure = plot(tse_all);
    title('TSE Complexity Against Trial')
    xlabel('Trial')
    ylabel('Complexity')
    grid on
    
    saveas(figure, (sprintf('/home/cs07/TSE_trial_se_%d.pdf',round(n))))
end
