%hvals = linspace(0, 1, 5);
%svals = linspace(0, 1, 5);

epsilonvals = [0.5, 0.7, 0.9]; %linspace(0, 1, 5);
Ccvals = [100, 500, 1000, 2000, 5000, 1e4]; %linspace(50, 1e3, 5);
hval = 0.3
sval = 0.95
nuvals = linspace(0, 0.1, 5)

%for hval = hvals
    for nuval = nuvals
        for epsilonval = epsilonvals
            for Ccval = Ccvals

                [newflog, P_time_D, P_time_W, newflogr, newPSD, newPSW]  = Full_PS_NHEJ(10, hval, sval, epsilonval, Ccval, nuval);

                home_directory = getenv('HOME');

                  round_h  = round(hval, 2);
                  round_s  = round(sval, 2);
                  round_ep = round(epsilonval, 2);                  
                  round_Cc = round(Ccval, 0);
round_nu = round(nuval, 2);

format_h  = sprintf('%.2f', round_h);
format_s  = sprintf('%.2f', round_s);
format_ep = sprintf('%.2f', round_ep);
format_Cc = sprintf('%.1f', round_Cc);
format_nu = sprintf('%.2f', round_nu);


                writematrix(newflog, sprintf("%s/FREQ_NHEJ/freq_time_h=%s_s=%s_ep=%s_Cc=%s_nu=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc, format_nu))
                writematrix(newflogr, sprintf("%s/FREQ/freq_space_h=%s_s=%s_ep=%s_Cc=%s_nu=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc, format_nu))
                
                writematrix(P_time_D, sprintf("%s/FFT_t_NHEJ/drive_h=%s_s=%s_ep=%s_Cc=%s_nu=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc, format_nu))
                writematrix(P_time_W, sprintf("%s/FFT_t_NHEJ/wild_h=%s_s=%s_ep=%s_Cc=%s_nu=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc, format_nu))
                
                writematrix(newPSD, sprintf("%s/FFT_k_NHEJ/drive_h=%s_s=%s_ep=%s_Cc=%s_nu=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc, format_nu))
                writematrix(newPSW, sprintf("%s/FFT_k_NHEJ/wild_h=%s_s=%s_ep=%s_Cc=%s_nu=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc, format_nu))
            end
        end
    end
%end

