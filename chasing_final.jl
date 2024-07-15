include("haplotype_model.jl")
using MAT
using Printf

hNval = 0.02
hval = 0.30
sval = 0.95

Cclist = [100, 500, 1000, 2000, 5000, 1e4]
eplist = [0.5, 0.7, 0.9]
nulist = range(0, 0.1, 5)
#nuval = 0
# for hval in range(0, stop=1, length=5)
#     println(hval)
     for nuval in nulist
        for epsilonval in eplist
            for Ccval in Cclist
                for i in 1:10

                    #sqmig_TSR (m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,width,height,D,pmig,preexisting,T,test_threshold,Gauss,terminate_resistance)
                    fem_freq, male_freq, population_size, ~, ~ = hexmig_TSR(1, sval, hval, hNval, 0.0, 6, Ccval*35^2, 0.0, 0.03, epsilonval, nuval, 0.0, 0.0, 0.0, 35, 35, 0.05, 0.0, 1000, 1000, 0, 0.0); #hexmig_TSR(1, 1, 0.3, 0.02, 0.01,6, 150*35^2,   0.0, 0.03, 0.995, 0.035, 0.0,0.0,0.0,35, 35, 0.05, 0.0, 300, 1000, 0, 0) # #hexmig_TSR(1, 1, 0.3, 0.02, 0.0, 6, 500, 0.0, 0.0, 0.995, 0.0, 0.0, 0.0, 0.0, 30, 30, 0.025, 0.0, 50, 1000, 0, 0.0);  
                    #hexmig_TSR(1, 1, 0.3, 0.02, 0.01,6, 150*35^2,   0.0, 0.03, 0.995, 0.0, 0.0,0.0,0.0,35, 35, 0.05, 0.0, 300, 1000, 0, 0)
                                                                            #(1, 1, 0.3, 0.02, 0.0, 6, carcap_fix, 0.0, ff, 0.95, 0.0, 0.0, 0.0, 0.0, valx, valy, diff_fix, pmig_fix, 0.0, TIME, 1000, 0, 0.0); #hexmig_TSR(1, 1, 0.3, 0.02, 0.0, 6, 500, 0.0, 0.0, 0.995, 0.0, 0.0, 0.0, 0.0, 30, 30, 0.025, 0.0, 50, 1000, 0, 0.0);  
                                                                            #(m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,width,height,pmig,preexisting,T,test_threshold,Gauss,terminate_resistance)
                    round_h  = round(hval, digits=2)
                    round_s  = round(sval, digits=2)
                    round_ep = round(epsilonval, digits=2)
                    round_Cc = round(Ccval, digits=0)
round_nu = round(nuval, digits=2)

format_h  = @sprintf("%.2f", round_h)
format_s  = @sprintf("%.2f", round_s)
format_ep = @sprintf("%.2f", round_ep)
format_Cc = @sprintf("%.1f", round_Cc)
format_nu = @sprintf("%.2f", round_nu)

                    home_directory = ENV["HOME"]
                    popmat_path_fem = joinpath(home_directory, "SIMS_NHEJ", "TM_fem_pop_$(i)_h=$(format_h)_s=$(format_s)_ep=$(format_ep)_Cc=$(format_Cc)_nu=$(format_nu).mat");
                    popmat_path_male = joinpath(home_directory,  "SIMS_NHEJ", "TM_male_pop_$(i)_h=$(format_h)_s=$(format_s)_ep=$(format_ep)_Cc=$(format_Cc)_nu=$(format_nu).mat"); #_nu=$(format_nu).mat")

                    population_size_mat  = float(population_size)
                    popmat_path = joinpath(home_directory, "SIMS_NHEJ", "TM_pop_$(i)_h=$(format_h)_s=$(format_s)_ep=$(format_ep)_Cc=$(format_Cc)_nu=$(format_nu).mat")


#println(population_size_mat[1, :, :, 1:5])
                    matwrite(popmat_path_male, Dict("population_male" => male_freq))
                    matwrite(popmat_path_fem, Dict("population_fem" => fem_freq))
                    matwrite(popmat_path, Dict("population" => population_size_mat))

                end
            end
        end
     end
# end


