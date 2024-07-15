hvals = linspace(0, 1, 5);
svals = linspace(0, 1, 5);
epsilonvals = linspace(0, 1, 5);
Ccvals = [100, 500, 1000, 2000, 5000, 1e4]; %linspace(50, 1e3, 5);
hval = 0.3;
sval = 0.95;

%for hval = hvals
%    for sval = svals
        for epsilonval = epsilonvals
            for Ccval = Ccvals

                [newflog, P_time_D, P_time_W, newflogr, newPSD, newPSW]  = Full_PS(10, hval, sval, epsilonval, Ccval);

P_time_W;
P_time_D;

                home_directory = getenv('HOME');

                  round_h  = round(hval, 2);
                  round_s  = round(sval, 2);
                  round_ep = round(epsilonval, 2);                  
                  round_Cc = round(Ccval, 0);

format_h  = sprintf('%.2f', round_h);
format_s  = sprintf('%.2f', round_s);
format_ep = sprintf('%.2f', round_ep);
format_Cc = sprintf('%.1f', round_Cc);



                writematrix(newflog, sprintf("%s/FREQ/freq_time_h=%s_s=%s_ep=%s_Cc=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc));
                writematrix(newflogr, sprintf("%s/FREQ/freq_space_h=%s_s=%s_ep=%s_Cc=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc));
                
                writematrix(P_time_D, sprintf("%s/FFT_t/drive_h=%s_s=%s_ep=%s_Cc=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc));
                writematrix(P_time_W, sprintf("%s/FFT_t/wild_h=%s_s=%s_ep=%s_Cc=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc));
                
                writematrix(newPSD, sprintf("%s/FFT_k/drive_h=%s_s=%s_ep=%s_Cc=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc));
                writematrix(newPSW, sprintf("%s/FFT_k/wild_h=%s_s=%s_ep=%s_Cc=%s.csv", home_directory, format_h, format_s, format_ep, format_Cc));
            end
        end
    %end
%end

