function [newflog, P_time_D, P_time_W, newflogr, newPSD, newPSW]  = Full_PS(num_files, hval, sval, epsilonval, Ccval)

format shortG

   FT_W = [];
   FT_D = [];
   nbin = 50; %number of blocks or points in PSD
   
   Fourier_space_W = zeros(num_files, nbin, nbin);
   Fourier_space_D = zeros(num_files, nbin, nbin);

                  round_h  = round(hval, 2);
                  round_s  = round(sval, 2);
                  round_ep = round(epsilonval, 2);                  
                  round_Cc = round(Ccval, 0);

format_h  = sprintf('%.2f', round_h);
format_s  = sprintf('%.2f', round_s);
format_ep = sprintf('%.2f', round_ep);
format_Cc = sprintf('%.1f', round_Cc)

    
    for i = 1:num_files
        home_directory = "/rds/general/user/csheeran/home"; %getenv('HOME');

        file_name1 = sprintf("%s/SIMS/TM_fem_pop_%d_h=%s_s=%s_ep=%s_Cc=%s.mat", home_directory, i, format_h, format_s, format_ep, format_Cc);
        load(file_name1);

        file_name2 = sprintf("%s/SIMS/TM_male_pop_%d_h=%s_s=%s_ep=%s_Cc=%s.mat", home_directory, i, format_h, format_s, format_ep, format_Cc);
        load(file_name2);
        
        file_name3 = sprintf("%s/SIMS/TM_pop_%d_h=%s_s=%s_ep=%s_Cc=%s.mat", home_directory, i, format_h, format_s, format_ep, format_Cc);
        load(file_name3);

        dim = size(population_male);
        Tall   = dim(4);
        time_chosen = round(3*Tall/4);
        
        T   = 701;
"POPULATION_FIN"
%population(1, :, :, 1:10)
    
        %TEMPORAL FT
    
        pop_W_av = reshape(mean(population_fem(:, :, 1, end-(T-1):end), [1, 2]), 1, []).*reshape(mean(population(1, :, :, end-(T-1):end), [2, 3]), 1, []) + reshape(mean(population_male(:, :, 1, end-(T-1):end), [1, 2]), 1, []).*reshape(mean(population(2, :, :, end-(T-1):end), [2, 3]), 1, [])
        pop_D_av = reshape(mean(population_fem(:, :, 4, end-(T-1):end), [1, 2]), 1, []).*reshape(mean(population(1, :, :, end-(T-1):end), [2, 3]), 1, []) + reshape(mean(population_male(:, :, 4, end-(T-1):end), [1, 2]), 1, []).*reshape(mean(population(2, :, :, end-(T-1):end), [2, 3]), 1, []);
        Fourier_time_W = fft(pop_W_av);
        Fourier_time_D = fft(pop_D_av);
    
        P_time_W  = Fourier_time_W .* conj(Fourier_time_W)/T;
        P_time_D  = Fourier_time_D .* conj(Fourier_time_D)/T;
    
        Ts = 1;          
        Fs = 1/Ts;               
        nyquist_frequency = Fs/2;               
    
        freq = (0:T-1)/T;
        mask = abs(freq) <= nyquist_frequency;
        freq = freq(mask);
        freq = freq(2:end); %remove DC
        flog = log(freq);

        P_time_W = 2*P_time_W(mask);
        P_time_W = P_time_W(2:end);
        P_time_D = 2*P_time_D(mask); %.* mask, times 2 as symmetric;
        P_time_D = P_time_D(2:end);
        
        [~, cent] = hist(flog,nbin); %use hist function to calc bin positions
        delta     = cent(2) - cent(1);
        
            for k = 1:nbin
                indp(k).indx  = find(flog >= cent(k)-delta/2 & flog < cent(k)+delta/2);
            end
        
        for k = 1:nbin    
            nblock(k)   = length(indp(k).indx);
            newflog(k)  = mean(freq(indp(k).indx));
            newPlogW(k) = mean(P_time_W(indp(k).indx));
            newPlogD(k) = mean(P_time_D(indp(k).indx));
            errlogW(k)  = newPlogW(k)/sqrt(nblock(k));
            errlogD(k)  = newPlogD(k)/sqrt(nblock(k));
        end

        nonnan_ind = find(~isnan(newflog));

        newflog  = newflog(nonnan_ind);
        newPlogW = newPlogW(nonnan_ind);
        newPlogD = newPlogD(nonnan_ind);

        FT_W = cat(1, FT_W, newPlogW);
        FT_D = cat(1, FT_D, newPlogD);
         

        
        %RADIAL
        [rows, cols] = size(population_male(:, :, 1, time_chosen));

        pop_W = population_fem(:, :, 1, time_chosen).*population(1, :, :, time_chosen) + population_male(:, :, 1, time_chosen).*population(2, :, :, time_chosen);
        pop_D = population_fem(:, :, 4, time_chosen).*population(1, :, :, time_chosen) + population_male(:, :, 4, time_chosen).*population(2, :, :, time_chosen);
        
        FT_space_W = fft2(pop_W);
        FT_space_D = fft2(pop_D);
        P_space_W  = FT_space_W .* conj(FT_space_W)/(rows*cols);
        P_space_D  = FT_space_D .* conj(FT_space_D)/(rows*cols);
        

        a       = 1; %2a of a hexagon
        a_prime = 2/(sqrt(3)*a);
        fs      = a_prime; 

        nyquist_frequency_row = fs/2;
        nyquist_frequency_col = fs/2;
        frequency_rows = (0:rows-1)/rows;
        frequency_cols = (0:cols-1)/cols;
        
        mask_rows = abs(frequency_rows) <= nyquist_frequency_row;
        mask_cols = abs(frequency_cols) <= nyquist_frequency_col;

        frequency_rows = frequency_rows(mask_rows);
        frequency_cols = frequency_cols(mask_cols);
        frequency_rows = frequency_rows(2:end);
        frequency_cols = frequency_cols(2:end);

        logrows = log(frequency_rows);
        logcols = log(frequency_cols);
        
        P_space_W = 2*P_space_W(mask_rows, mask_cols); 
        P_space_D = 2*P_space_D(mask_rows, mask_cols); 
        P_space_W = P_space_W(2:end, 2:end);
        P_space_D = P_space_D(2:end, 2:end);
        
        [~, centr] = hist(logrows, nbin); %use hist function to calc bin positions
        deltar     = centr(2) - centr(1);

        [~, centc] = hist(transpose(logcols), nbin); %use hist function to calc bin positions
        deltac     = centc(2) - centc(1);
        
        
            for k = 1:nbin
                indpx(k).indx  = find(logrows >= centr(k)-deltar/2 & logrows < centr(k)+deltar/2);
                indpy(k).indx  = find(logcols >= centc(k)-deltac/2 & logcols < centc(k)+deltac/2);
            end

        for k = 1:nbin
            newflogx(k)  = mean(frequency_rows(indpx(k).indx));
            for l = 1:nbin
                newflogy(l)  = mean(frequency_cols(indpy(l).indx));
            
                meanW = nanmean(P_space_W(indpx(k).indx, indpy(l).indx));
                meanD = nanmean(P_space_D(indpx(k).indx, indpy(l).indx));

            if size(meanW, 2) > 1
                meanW = nanmean(meanW);
            end

            if size(meanD, 2) > 1
                meanD = nanmean(meanD);
            end

            if size(meanW, 2) < 1
                meanW = nan;
            end

            if size(meanD, 2) < 1
                meanD = nan;
            end

            newPlogW(k, l) = meanW;
            newPlogD(k, l) = meanD;
            end
        end
       

        Fourier_space_W(i, :, :) = newPlogW;
        Fourier_space_D(i, :, :) = newPlogD;


    end

    %TEMPORAL:
    P_time_W = mean(FT_W, 1);
    P_time_D = mean(FT_D, 1);

    %newflog
    
    %RADIAL:
    nonnan_indx = find(~isnan(newflogx));
    nonnan_indy = find(~isnan(newflogy));
    newflogx   = newflogx(nonnan_indx);
    newflogy   = newflogy(nonnan_indy);

    Fourier_space_D = Fourier_space_D(:, nonnan_indx, nonnan_indy);
    Fourier_space_W = Fourier_space_W(:, nonnan_indx, nonnan_indy);

    P_space_D = mean(Fourier_space_D, 1);
    P_space_D = reshape(P_space_D, [size(Fourier_space_D, 2), size(Fourier_space_D, 3)]);

    P_space_W = mean(Fourier_space_W, 1);
    P_space_W = reshape(P_space_W, [size(Fourier_space_W, 2), size(Fourier_space_W, 3)]);
   
    ind = 0;
    for i=1:size(newflogx, 2)
        for j=1:size(newflogy, 2)
            ind = ind +1;
            newflogr(ind) = sqrt(newflogx(i)^2 + newflogy(j)^2);
            newPSD(ind) = P_space_D(i, j);
            newPSW(ind) = P_space_W(i, j);
        end
    end

%newflogr, newPSD, newPSW


end

