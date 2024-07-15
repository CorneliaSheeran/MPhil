using Distributions
using LinearAlgebra
using PoissonRandom #probably not needed for functions in this file
using StatsBase
using SpecialFunctions

function linindices(m,n)
    # Returns the ((n-1)^m+1)xm array [1 1 1 ... 1;2 1 ... 1;3 1 ... 1; ... ; n-1 1 ... 1;2 1 ... 1;2 2 1 ... 1; ... ;n-1 n-1 ... n-1; n n n]
    # Try, say, linindices(3,4) 

    na = (n-1)^m+1;

    linindicesarray = zeros(Int,na,m);

    for i = 1:na-1
        for j = 1:m
            linindicesarrayindex = ceil(Int, i/((n-1)^(j-1)));
            linindicesarray[i,j] = mod(linindicesarrayindex, (n-1));

            if linindicesarray[i,j] == 0
                linindicesarray[i,j] = n-1;
            end

        end
    end

    for j = 1:m
        linindicesarray[na,j] = n;
    end

    return linindicesarray
end

function femalefitnessmatrix(m,n,s,hD,hN,σ)
    # Outputs the matrix Wf, where Wf[i,j] is the fitness of a female with genotype i,j

    na = (n-1)^m + 1;
    #print("na = ")
    #print(na)
    W = zeros(na);
    Wf = zeros(na,na);

    indexarray = linindices(m,n);
    #print("indexarray = ")
    #print(indexarray)

    nW = zeros(Int, na-1);
    nR = zeros(Int, na-1); #or na??

    for k = 1:na-1
        #print("haplotype = ")
        haplotypevector = indexarray[k,:];
        #print(haplotypevector)

        nn = zeros(Int, n);#changed from n-1

        for kk = 1:m
            nn[haplotypevector[kk]] +=1;
        end
        #print("number of alleles = ")
        #print(nn)
        if nn[3] == 0
            nW[k] = nn[1];
            #print("nW[k] = ")
            #print(nW[k])
            nR[k] = nn[2];
            #print("nR[k] = ")
            #print(nR[k])
        else nR[k] = -1;

        end

        w = 1;

        if sum(nn[3:end]) == 0 #no deleterious alleles present
            for kk = 1:m
                if haplotypevector[kk] == 2
                    w = w * (1 - σ[kk]);
                end
            end

        else w = 1-s;
        end

        W[k] = w;
        #print("W[k] = ")
        #print(W[k])
    end

    W[na] = 1-s;

    for i = 1:na
        for j = 1:na
            #print("i = ")
            #print(i)
            #print("j = ")
            #print(j)
            deleterious = zeros(Int, 2);
            if i == na || nR[i] == -1
                #print("yes for i")
                deleterious[1] = 1;
            end
            if j == na || nR[j] == -1
                #print("yes for j")
                deleterious[2] += 1;
                #print("deleterious[2] = ")
                #print(deleterious[2])
            end
            #print("deleterious?")
            #print(deleterious)

    # QUESTION: Do the effects of multiple N alleles across target sites compound?

            if deleterious[1] == 1 && deleterious[2] == 1
                Wf[i,j] = 1-s;

            elseif i == na
                #before, this was nR[j], without "m - ". It had hD rather than hN.
                #The idea was that a single R allele confers partial resistance.
                #Now, all have to be R.
                if m - nR[j] == 0 #First heterozygous R...R/D...D
                    Wf[i,j] = (1 - hN*s) * W[j];
                elseif m - nR[j] > 0 #If any W alleles present, full cost
                                     #imposed
                    Wf[i,j] = (1 - hD*s) * W[j];
                end

            elseif j == na

                if m - nR[i] == 0
                    Wf[i,j] = (1 - hN*s) * W[i];
                elseif m - nR[i] > 0
                    Wf[i,j] = (1 - hD*s) * W[i]
                end

            elseif nR[i] == -1
                Wf[i,j] = (1 - hN*s) * W[j];
            elseif nR[j] == -1
                Wf[i,j] = (1 - hN*s) * W[i];
            else Wf[i,j] = W[i]*W[j];
            end

        end

    end

    return Wf

end

function drivematrix(m,n,ϵ,ν,β)
    # Prints out the matrix κ. (I+κ)[i,j] is the probability of transition j->i under drive (ie excluding de novo mutations).

    na = (n-1)^m + 1;

    #kappasingle[i,j] is the probability of transition j/D->i/D during gametogenesis for a j/D heterozygote
    kappasingle = [1-ϵ 0.0 0.0 0.0;ϵ*ν*β 1.0 0.0 0.0;ϵ*ν*(1-β) 0.0 1.0 0.0;ϵ*(1-ν) 0.0 0.0 1.0];

    #κ[i,j] is the probability of transition j->i during gametogenesis (multilpex)
    κ = zeros(na,na)

    indexarray = linindices(m,n);

    for j=1:na-1

        #Introduce variable to keep track of total column probability, for use
        #later
        probtracker = 0.0;
        for i=1:na-1
            for kk=1:m
                #Calculate alleles at site m=kk for haplotypes i,j
                #print("ii = ")
                ii = indexarray[i,kk];
                #print(ii)
                #print("jj = ")
                jj = indexarray[j,kk];
                #print(jj)

                #First, the m=1 target site allele must transition j->i
                if kk==1
                    #print("kk = 1")
                    κ[i,j] = kappasingle[ii,jj];
                    #print("κ")
                    #print(i)
                    #print(j)
                    #print("=")
                    #print(κ[i,j])
                #Next, each subsequent allele must transition as well
                else
                    #print("kk =")
                    #print(kk)
                    κ[i,j] = κ[i,j]*kappasingle[ii,jj];
                    #print("κ")
                    #print(i)
                    #print(j)
                    #print("=")
                    #print(κ[i,j])
                end

            end

            probtracker = probtracker + κ[i,j];

        end

        #Finally, P(j->drive)=1-sum_{i<drive}P(j->i)
        κ[na,j] = 1 - probtracker;
    end

    #P(drive -> drive) = 1
    κ[na,na] = 1.0;

    #Subtract identity matrix
    for i=1:na
        κ[i,i] = κ[i,i] - 1;
    end

    return κ

end

function mutationmatrix(m,n,μ,ξ)
    #Outputs mutmatrix. mutmatrix[i,j] denotes the probability of de novo mutation from
    #haplotype j to haplotype i.

    #First calculate the number of (ordered) haplotypes
    na = (n-1)^m + 1;

    #muvec = [1 - muSite; xi * muSite; (1 - xi) * muSite];

    #uR = [0; 1; 0];
    #uN = [0; 0; 1];

    #singlesitemut is the matrix of mutation probabilities j->i for m=1
    singlesitemut = [1-μ 0.0 0.0 0.0;ξ*μ 1.0 0.0 0.0;(1-ξ)*μ 0.0 1.0 0.0;0.0 0.0 0.0 1.0];

    #Probability of mutation from j to i
    mutmatrix = zeros(na,na);

    indexarray = linindices(m,n);

    for i = 1:na-1
        for j = 1:na-1
            for kk = 1:m
                ii = indexarray[i,kk];
                jj = indexarray[j,kk];
                #all entries are initially zero, so have to set first value
                if kk == 1
                    mutmatrix[i,j] = singlesitemut[ii,jj];
                else
                    mutmatrix[i,j] = mutmatrix[i,j]*singlesitemut[ii,jj];
                end
            end
        end
    end

    mutmatrix[na,na] = 1.0;

    #subtract identity matrix??

    return mutmatrix;

end

function newfrequency(Rm,α,xold,yold,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss)
    # Input Rm (growth rate), α (function of Rm and carrying capacity)
    # xold (vector of old female allele frequencies)
    # yold (vector of old male allele frequencies)
    # μ not actually needed
    # Wf female fitness matrix, from femalefitnessmatrix.jl
    # Wm male fitness matrix - for our simulations this is always a matrix of 1s
    # na (number of haplotypes)
    # ψ (probability that an individual offspring is male)
    # mutmatrix probability matrix of de novo mutations, from mutationmatrix.jl
    # κ matrix of drive induced transitions, from drivematrix.jl
    # Nf old female population size
    # Nm old male population size
    # Gauss = 1 implements Gauss-Poisson hybrid approximation of multinomial distribution. Gauss = 0 uses native Julia program.

    # Outputs:
    # (xnew,ynew,Nfnew,Nmalenew,Nnew) = (new female haplotype frequency vector, new male haplotype frequency vector, new female population, new male population, new total population, ie Nf + Nmale)

    #if either sex has died out, there will be no reproduction
    if (mean(xold) == 0 || mean(yold) == 0) || (Nf == 0 || Nm == 0)
        xnew = zeros(na);
        ynew = zeros(na);
        return (xnew,ynew,0,0,0)
    end
    #First calculate average fitnesses
    meanf = transpose(xold)*Wf*yold;
    meanm = transpose(xold)*Wm*yold;
    if meanf == 0 || meanm == 0
        xnew = zeros(na);
        ynew = zeros(na);
        return (xnew,ynew,0,0,0)
    end

    #Now calculate the new frequencies obtained by random mating and ne novo
    #mutation. At this stage we NEGLECT DRIVE.
    xnewnodrive = mutmatrix * (((Wf*yold).*xold + (Wf*xold).*yold)/(2*meanf));
    ynewnodrive = mutmatrix * (((Wm*xold).*yold + (Wm*yold).*xold)/(2*meanm));

    #Now correct for drive
    xdriveterm = zeros(na);
    ydriveterm = zeros(na);

    for i = 1:na
        for j = 1:na
            #Check it's OK to allow all haplotypes. #Yes, I think it is.
            #H double counts drive homozygotes, but this is OK since κ zeros
            #them anyway
            #This has been modified from paper. Multiply by mutmatrix[j,j]
            #instead of (1-μ)
            H = xold[j]*yold[na] + xold[na]*yold[j];
            xdriveterm[i] += 0.5*Wf[j,na]*mutmatrix[j,j]*κ[i,j]*H/meanf;
            ydriveterm[i] += 0.5*Wm[j,na]*mutmatrix[j,j]κ[i,j]*H/meanm;
        end
        #print("xdriveterm = ")
        #println(xdriveterm)
        #print("ydriveterm = ")
        #println(ydriveterm)
    end

    xoverall = xnewnodrive + xdriveterm;
    if minimum(xoverall) < 0
        print("error - new expected female frequency has negative entry")
        return (-1,-1,-1,-1,-1)
    end
    yoverall = ynewnodrive + ydriveterm;
    if minimum(yoverall) < 0
        print("error - new expected male frequency has negative entry")
        return (-1,-1,-1,-1,-1)
    end
    #print("xoverall ")
    #println(xoverall)
    #print("yoverall ")
    #println(yoverall)

    #Calculate expected population size in next generation according to
    #Beverton-Holt
    Nexptd = 2*Rm*meanf*Nf/(1 + (Nf+Nm)/α);
    if Nexptd < 0 || isnan(Nexptd) == 1
        print("Nexptd = ")
        print(Nexptd)
        print("mean female fitness = ")
        print(meanf)
        print("Number of females = ")
        print(Nf)
        print("Number of males = ")
        print(Nm)
        print("α = ")
        print(α)
    end
    #Nnew = rand(Poisson(Nexptd));

    #Now sample according to Wright-Fisher
    if Nexptd > 10^15
        #Check with supervisor!!!
        Nmalenew = round(Int128,rand(Normal(ψ*Nexptd,sqrt(ψ*Nexptd))));
        Nfnew = round(Int128,rand(Normal((1-ψ)*Nexptd,sqrt((1-ψ)*Nexptd))));
    else
        Nmalenew = pois_rand(ψ*Nexptd);
        Nfnew = pois_rand((1-ψ)*Nexptd);
    end
    Nnew = Nmalenew + Nfnew;
    #Nmalenew = rand(Binomial(Nnew,ψ));
    #Nfnew = Nnew - Nmalenew;
    ###Alternatively, Nmalenew~Poi(ψ*Nexptd) and Nfnew~Poi((1-ψ)*Nexptd);
    if Gauss == 0
        xnumnew = rand(Multinomial(Nfnew,xoverall));
        ynumnew = rand(Multinomial(Nmalenew,yoverall));
    else
        xnumnew = GaussPoissonHybrid_mnrnd(Nfnew,xoverall);
        ynumnew = GaussPoissonHybrid_mnrnd(Nmalenew,yoverall);
    end

    if Nfnew > 0
        xnew = xnumnew/Nfnew;
    else
        xnew = zeros(na);
    end

    if Nmalenew > 0
        ynew = ynumnew/Nmalenew;
    else
        ynew = zeros(na);
    end

    return (xnew,ynew,Nfnew,Nmalenew,Nnew)
end

function resistancehaplotypes(m,n,resistancealleles)
    # if resistancealleles is a vector of the resistant allele positions,
    # resistancehaplotypes calculates the resistant haplotype positions
        # eg. if R,N are resistant and m=2, RR,RN,NR,NN are resistant haplotypes
    
        na = (n-1)^m + 1;
        resistancehaplotypes = zeros(Int,na);
    
        indexarray = linindices(m,n);
    
    
        for iiii = 1:na
            checker = 1;
            for jjjj = 1:m
                if (indexarray[iiii,jjjj] in resistancealleles) == false
                    checker = 0;
                end
            end
    
            if checker == 1
                resistancehaplotypes[iiii] = 1;
            end
        end
    
        return resistancehaplotypes
    
end

function TargetSiteResistance(m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,preexisting,T,test_threshold,Gauss)
    # m = number of target sites
    # s = homozygous fitness cost of drive
    # s*h = fitness cost of heterozygous W/D individual
    # s*hN = fitness cost of heterozygous W/N or R/D individual
    # σ_i = fitness cost of R allele at target site i, 1≦i≦m
    # Rm =
    # K =
    # xD,yD = initial proportion of D...D haplotypes in females/males respectively
    # ϵ = efficiency of drive cleavage per target site
    # ν = non-homologous end-joining (NHEJ) rate per generation per individual per
    # target site
    # μ = mutation rate per cleavage target site per generation per individual per
    # target site (e.g. if a target site contains 10bp and bp mutation rate is 1e-9,
    # μ=1e-8)
    # β = fraction of functional resistant NHEJ mutations
    # ξ = fraction of functional resistant single nucleotide mutations
    # preexisting = 1,0 to run a burn in period of 1/σ and set initial allele
    # frequencies or not, respectively
    # T = length of simulations in generations
    
        #Set proportion of males
        ψ = 0.5;
    
        N=round(Int128,K);#Make this better later
        if Gauss == 0 && log2(N) > 63
            print("Error - if population > 2^63 then must use Gauss-Poisson Approximation")
            return 2
        end
        #Nm=ψ*N;
        #Nf=N-Nm;
        if Gauss == 0
            Nm = rand(Binomial(N,ψ));
            Nf = N - Nm;
        else
            vvv = GaussPoissonHybrid_mnrnd(N,[ψ;1-ψ]);
            Nm = round(Int128,vvv[1]);
            Nf = round(Int128,vvv[2]);
        end
        #Density dependent parameter from Beverton-Holt model of population dynamics
        α = K/(Rm - 1);
    
        #Dominance coefficients
        hD = h; #Why not just input hD instead of h and avoid this???
    
        if length(σ) == 1
            σ = σ*fill(1.0,(m,1));
        elseif length(σ) != m
            print("Error - σ should be scalar or vector of length m")
            return
        end
    
        #Set number of alleles
        n = 4;
    
        #Calculate number of (ordered) haplotypes
        na = (n-1)^m + 1;
    
        #Calculate fitness matrices
        Wm = fill(1.0,(na,na));
        Wf = femalefitnessmatrix(m,n,s,hD,hN,σ);
    
        #Calculate mutation probaility matrix. M[i,j] is the probability of de novo
        #mutation haplotype j -> haplotype i
        mutmatrix = mutationmatrix(m,n,μ,ξ);
    
        #Calculate drive conversion matrix. Suppose an individual has genotype
        #[j,drive]. Then K[i,j]+δ_{ij} is the probability of the conversion j->i via
        #drive and NHEJ mutations. So K[i,j] is the mean relative change in
        #frequency j->i.
        κ = drivematrix(m,n,ϵ,ν,β);
    
        ## Initial frequencies
        x0 = zeros(na);
        x0[1] = 1.0;
        y0 = zeros(na);
        y0[1] = 1.0;
    
        ##If preexisting
        if preexisting == 1 && μ > 0
            #Calculate burn in time
            if minimum(σ) <= 0
                print("Error - for preexisting = 1, resistance alleles must always incur fitness cost")
                return
            end
    
            Tburn = round(1/minimum(σ));
    
            # Now initialise as in original code
            indm = mutmatrix.>0;
            for i = 1:na
                for j = 2:na
                    indm[i,j] = 0;
                end
            end
            indm[1,1] = 0;
    
            muin = mutmatrix[indm];
            #print("muin = ")
            #println(muin)
    
            ffff = diagm((Wf[indm].-1)./2);
            #print("ffff = ")
            #println(ffff)
            #println(ffff)
    
            indm2 = mutmatrix.>0;
            indm2[1,1] = 0;
            for i = 2:na
                for j = 2:na
                    if indm2[i,1] > 0 && indm2[j,1] > 0
                        indm2[i,j] = 1;
                    else
                        indm2[i,j] = 0;
                    end
                end
            end
    
            for i = 1:na
                indm2[i,1] = 0;
                indm2[1,i] = 0;
            end
    
            modmut1 = mutmatrix[indm2];
            lengthvec = length(modmut1);
            if lengthvec != length(ffff);
                println("Error with preexisting")
                println(lengthvec)
                println(length(ffff))
            end
            lengthvec = round(Int,sqrt(lengthvec));
            modmut2 = zeros((lengthvec,lengthvec));
            for i = 1:lengthvec
                for j = 1:lengthvec
                    modmut2[i,j] = modmut1[lengthvec*(j-1)+i];
                    if i == j
                        modmut2[i,j] += -1;
                    end
                end
            end
            #print("modmut2 = ")
            #println(modmut2)
    
            indmvec = BitArray(undef,na);
            for i = 1:na
                indmvec[i] = indm[i,1];
            end
    
            if det(modmut2+ffff) != 0
            x0[indmvec] = -inv(modmut2+ffff)*muin;
            y0[indmvec] = -inv(modmut2+ffff)*muin;
            else
                Tburn = Tburn*10;
            end
            #println(x0)
            #println(y0)
    
            x0[na] = 0;
            y0[na] = 0;
            x0[1] = 1 - sum(x0[2:end]);
            y0[1] = 1 - sum(y0[2:end]);
            #print("x0 = ")
            #println(x0)
    
            for tt = 1:Tburn
                #m
                (x0,y0,Nf,Nm,N) = newfrequency(Rm,α,x0,y0,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss);
                if (x0,y0,Nf,Nm,N) == (-1,-1,-1,-1,-1)
                    println("error occured during burn in, time = ")
                    print(tt)
                    return
                end
            end
            #print(x0)
            #print(y0)
            #print(N)
        end
    
        ## Now run simulation, keeping track of frequencies and population size
        # First, add drive to initial frequencies
        x0[na] = xD;
        x0[1] = x0[1] - xD;
        y0[na] = yD;
        y0[1] = y0[1] - yD;
    
        if x0[1] < 0 || y0[1] < 0
            print("error - not enough wild alleles to add drive")
            return 5
        end
    
        xcurrent = x0;
        ycurrent = y0;
        # x,y will keep track over time of fe/male haplotype frequencies,
        # respectively
        x = zeros((na,T));
        y = zeros((na,T));
        time = zeros(Int,T);
        x[:,1] = xcurrent;
        y[:,1] = ycurrent;
    
        #popsize keeps track over time of female, male and total populations
        popsize = zeros(Int128, (3,T));
        popsize[:,1] = [Nf;Nm;N];
    
        #We want to keep track of the resistant allele
        testallele = 2;
        resistance = 0;
        resistancepositions = resistancehaplotypes(m,n,[2,3]);
        #Calculate the number of resistance alleles in each haplotype
        indexarray = linindices(m,n);
    
        for t = 1:T-1
            time[t+1] = t;
            (xcurrent,ycurrent,Nf,Nm,N) = newfrequency(Rm,α,xcurrent,ycurrent,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss);
            if (xcurrent,ycurrent,Nf,Nm,N) == (-1,-1,-1,-1,-1)
                print("error occured during main phase, time = ")
                print(t)
                return(x,y,popsize)
            end
            x[:,t+1] = xcurrent;
            y[:,t+1] = ycurrent;
            popsize[:,t+1] = [Nf;Nm;N];
    
            # Test for resistance
    
            resistfreqf = dot(xcurrent,resistancepositions);
            resistfreqm = dot(ycurrent,resistancepositions);
            #resistancefreq = (1/N)*(Nf*resistfreqf + Nm*resistfreqm);
            if resistfreqf >= test_threshold && resistfreqm >= test_threshold
                if resistance == 0
                    τ = t;
                    #print(xcurrent)
                    #print(ycurrent)
                    #print(resistfreqf)
                    #print(resistfreqm)
                    #print(resistancepositions)
                end
                resistance = 1;
            end
    
        end
    
        if resistance == 0
            τ = NaN;
        end
    
        #recovery = 0;
        #if popsize[3,T] >= 0.1*K
        #    recovery = 1;
        #end
    
        return (x,y,popsize,resistance,τ)
end

function square_deme_migration(A,pmig,Gauss)
    # The entry A[i,j] denotes the initial population of individuals in deme [i,j].
    # Individuals migrate between dquare demes
    # Each individual migrates with probability pmig. The destination is chosen
    # uniformly randomly from the neighbouring demes
        acrm = size(A)[1];
        upn = size(A)[2];
    
        # This program only works for at least 2x2 matrices
        if acrm<2 || upn<2
            print("Error - need at least 2x2 grid for square demes")
            return 1
        end
        #println("α")
    
        # Create arrays of the number of individuals emigrating from each deme.
        # For example, remain[i,j] denotes the number of individuals moving
        # [i,j]->[i,j], whereas down[i,j]  represents the number moving
        # [i,j]->[i,j-1].
        if Gauss == 0
            remain = zeros(Int,(acrm,upn));
            up = zeros(Int,(acrm,upn));
            down = zeros(Int,(acrm,upn));
            left = zeros(Int,(acrm,upn));
            right = zeros(Int,(acrm,upn));
        else
            remain = zeros(Int128,(acrm,upn));
            up = zeros(Int128,(acrm,upn));
            down = zeros(Int128,(acrm,upn));
            left = zeros(Int128,(acrm,upn));
            right = zeros(Int128,(acrm,upn));
        end
        #println("β")
    
        for i = 1:acrm
            for j = 1:upn
                if Gauss == 0
                    # Need population to have type Int
                    A_temp = zeros(Int,(acrm,upn));
                    for i = 1:acrm
                        for j = 1:upn
                            A_temp[i,j] = trunc(Int,A[i,j]);
                        end
                    end
                    A = A_temp;
                    # First consider interior demes
                    if 1<i<acrm && 1<j<upn
                        #println("γ")
                        v = rand(Multinomial(A[i,j],[1-pmig;pmig/4;pmig/4;pmig/4;pmig/4]));
                        #println("δ")
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        up[i,j] = v[3];
                        left[i,j] = v[4];
                        down[i,j] = v[5];
    
                    # Now non-corner demes on the left edge
                    elseif i == 1 && 1<j<upn
                        v = rand(Multinomial(A[i,j],[1-(3/4)*pmig;pmig/4;pmig/4;pmig/4]));
    
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        up[i,j] = v[3];
                        # no left movement
                        down[i,j] = v[4];
    
                    # Now right edge
                    elseif i == acrm && 1<j<upn
                        v = rand(Multinomial(A[i,j],[1-(3/4)*pmig;pmig/4;pmig/4;pmig/4]));
    
                        remain[i,j] = v[1];
                        up[i,j] = v[2];
                        left[i,j] = v[3];
                        down[i,j] = v[4];
    
                    # Now bottom edge
                    elseif 1<i<acrm && j == 1
                        v = rand(Multinomial(A[i,j],[1-(3/4)*pmig;pmig/4;pmig/4;pmig/4]));
    
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        up[i,j] = v[3];
                        left[i,j] = v[4];
    
                    elseif 1<i<acrm && j == upn
                        v = rand(Multinomial(A[i,j],[1-(3/4)*pmig;pmig/4;pmig/4;pmig/4]));
    
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        left[i,j] = v[3];
                        down[i,j] = v[4];
    
                    elseif i == 1 && j == 1
                        #println("γ")
                        v = rand(Multinomial(A[i,j],[1-(1/2)*pmig;pmig/4;pmig/4]));
                        #println("δ")
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        up[i,j] = v[3];
    
                    elseif i == acrm && j == 1
                        v = rand(Multinomial(A[i,j],[1-(1/2)*pmig;pmig/4;pmig/4]));
                        remain[i,j] = v[1];
                        up[i,j] = v[2];
                        left[i,j] = v[3];
    
                    elseif i == 1 && j == upn
                        v = rand(Multinomial(A[i,j],[1-(1/2)*pmig;pmig/4;pmig/4]));
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        down[i,j] = v[3];
    
                    elseif i == acrm && j == upn
                        v = rand(Multinomial(A[i,j],[1-(1/2)pmig;pmig/4;pmig/4]));
                        remain[i,j] = v[1];
                        left[i,j] = v[2];
                        down[i,j] = v[3];
                    else
                        print("Error calculating emigrant numbers")
                        return 2
                    end
                else
                    # Need population to have type Int128
                    A_temp = zeros(Int128,(acrm,upn));
                    for i = 1:acrm
                        for j = 1:upn
                            A_temp[i,j] = trunc(Int,A[i,j]);
                        end
                    end
                    A = A_temp;
                    # First consider interior demes
                    if 1<i<acrm && 1<j<upn
                        #println("γ")
                        v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/4;pmig/4;pmig/4;pmig/4]);
                        #println("δ")
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        up[i,j] = v[3];
                        left[i,j] = v[4];
                        down[i,j] = v[5];
    
                    # Now non-corner demes on the left edge
                    elseif i == 1 && 1<j<upn
                        v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/3;pmig/3;pmig/3]);
    
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        up[i,j] = v[3];
                        # no left movement
                        down[i,j] = v[4];
    
                    # Now right edge
                    elseif i == acrm && 1<j<upn
                        v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/3;pmig/3;pmig/3]);
    
                        remain[i,j] = v[1];
                        up[i,j] = v[2];
                        left[i,j] = v[3];
                        down[i,j] = v[4];
    
                    # Now bottom edge
                    elseif 1<i<acrm && j == 1
                        v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/3;pmig/3;pmig/3]);
    
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        up[i,j] = v[3];
                        left[i,j] = v[4];
    
                    elseif 1<i<acrm && j == upn
                        v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/3;pmig/3;pmig/3]);
    
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        left[i,j] = v[3];
                        down[i,j] = v[4];
    
                    elseif i == 1 && j == 1
                        #println("γ")
                        v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/2;pmig/2]);
                        #println(v)
                        #println("δ")
                        remain[i,j] = v[1];
                        #println("δ.5")
                        right[i,j] = v[2];
                        up[i,j] = v[3];
                        #println("ϵ")
    
                    elseif i == acrm && j == 1
                        v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/2;pmig/2]);
                        remain[i,j] = v[1];
                        up[i,j] = v[2];
                        left[i,j] = v[3];
    
                    elseif i == 1 && j == upn
                        v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/2;pmig/2]);
                        remain[i,j] = v[1];
                        right[i,j] = v[2];
                        down[i,j] = v[3];
    
                    elseif i == acrm && j == upn
                        v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/2;pmig/2]);
                        remain[i,j] = v[1];
                        left[i,j] = v[2];
                        down[i,j] = v[3];
                    else
                        print("Error calculating emigrant numbers")
                        return 2
                    end
                end
            end
        end
    
        for i = 1:acrm
            if down[i,1] != 0
                print("Error calculating emigrant numbers")
                return 3
            elseif up[i,upn] != 0
                print("Error calculating emigrant numbers")
                return 4
            end
        end
    
        for j = 1:upn
            if left[1,j] != 0
                print("Error calculating emigrant numbers")
                return 5
            elseif right[acrm,j] != 0
                print("Error calculating emigrant numbers")
                return 6
            end
        end
    
        # Now calculate population distribution after migration, B[i,j]
        if Gauss == 0
            B = zeros(Int,(acrm,upn));
        else
            B = zeros(Int128,(acrm,upn));
        end
    
        for i = 1:acrm
            for j = 1:upn
                if 1<i<acrm && 1<j<upn
                    B[i,j] = remain[i,j] + right[i-1,j] + up[i,j-1] + left[i+1,j] + down[i,j+1];
                elseif i == 1 && 1<j<upn
                    B[i,j] = remain[i,j] + up[i,j-1] + left[i+1,j] + down[i,j+1];
                elseif i == acrm && 1<j<upn
                    B[i,j] = remain[i,j] + right[i-1,j] + up[i,j-1] + down[i,j+1];
                elseif 1<i<acrm && j == 1
                    B[i,j] = remain[i,j] + right[i-1,j] + left[i+1,j] + down[i,j+1];
                elseif 1<i<acrm && j == upn
                    B[i,j] = remain[i,j] + right[i-1,j] + up[i,j-1] + left[i+1,j];
                elseif i == 1 && j == 1
                    B[i,j] = remain[i,j] + left[i+1,j] + down[i,j+1];
                elseif i == acrm && j == 1
                    B[i,j] = remain[i,j] + right[i-1,j] + down[i,j+1];
                elseif i == 1 && j == upn
                    B[i,j] = remain[i,j] + up[i,j-1] + left[i+1,j];
                elseif i == acrm && j == upn
                    B[i,j] = remain[i,j] + right[i-1,j] + up[i,j-1];
                end
            end
        end
    
        return B
    
end

function hex_deme_migration(A,pmig)
    # The entry A[i,j] denotes the initial population of individuals in deme [i,j].
    # Individuals migrate between dquare demes
    # Each individual migrates with probability pmig. The destination is chosen
    # uniformly randomly from the neighbouring demes
        acrm = size(A)[1];
        upn = size(A)[2];
    
        # This program only works for at least 2x2 matrices
        if acrm<2 || upn<2
            print("Error - need at least 2x2 grid for square demes")
            return 1
        end
        #println("α")
    
        # Create arrays of the number of individuals emigrating from each deme.
        # For example, remain[i,j] denotes the number of individuals moving
        # [i,j]->[i,j], whereas down[i,j]  represents the number moving
        # [i,j]->[i,j-1].

        remain = zeros(Int,(acrm,upn));
        up = zeros(Int,(acrm,upn));
        down = zeros(Int,(acrm,upn));
        southwest = zeros(Int,(acrm,upn));
        southeast = zeros(Int,(acrm,upn));
        northwest = zeros(Int,(acrm,upn));
        northeast = zeros(Int,(acrm,upn));

        #println("β")

        # Need population to have type Int
        A_temp = zeros(Int,(acrm,upn));
        for i = 1:acrm
            for j = 1:upn
                A_temp[i,j] = trunc(Int,A[i,j]);
            end
        end
        A = A_temp;
    
        for i = 1:acrm
            for j = 1:upn
                    
                    # First consider interior demes
                    if 1<i<acrm && 1<j<upn
                        #println("γ")
                        v = rand(Multinomial(A[i,j],[1-pmig;pmig/6;pmig/6;pmig/6;pmig/6;pmig/6;pmig/6]));
                        #println("δ")
                        remain[i,j] = v[1];
                        up[i,j] = v[2];
                        down[i,j] = v[3];
                        southwest[i,j] = v[4];
                        southeast[i,j] = v[5];
                        northwest[i,j] = v[6];
                        northeast[i,j] = v[7];
    
                    # Now non-corner demes on the left edge
                    elseif i == 1 && 1<j<upn
                        v = rand(Multinomial(A[i,j],[1-(4/6)*pmig;pmig/6;pmig/6;pmig/6;pmig/6]));
    
                        remain[i,j] = v[1];
                        up[i,j] = v[2];
                        down[i,j] = v[3];
                        # no southwest movement
                        southeast[i,j] = v[4];
                        # no northwest movement
                        northeast[i,j] = v[5];
    
                    # Now right edge
                    elseif i == acrm && 1<j<upn
                        v = rand(Multinomial(A[i,j],[1-(4/6)*pmig;pmig/6;pmig/6;pmig/6;pmig/6]));
    
                        remain[i,j] = v[1];
                        up[i,j] = v[2];
                        down[i,j] = v[3];
                        southwest[i,j] = v[4];
                        # no southeast movement
                        northwest[i,j] = v[5];
                        # no northeast movement
    
                    # Now bottom edge
                    elseif 1<i<acrm && j == 1
                        if i%2 == 1 # Odd i index hexagons are different to even i index hexagons
                            v = rand(Multinomial(A[i,j],[1-(3/6)*pmig;pmig/6;pmig/6;pmig/6]));

                            remain[i,j] = v[1];
                            up[i,j] = v[2];
                            # no down movement
                            # no southwest movement
                            # no southeast movement
                            northwest[i,j] = v[3];
                            northeast[i,j] = v[4];

                        else
                            v = rand(Multinomial(A[i,j],[1-(5/6)*pmig;pmig/6;pmig/6;pmig/6;pmig/6;pmig/6]));
                            remain[i,j] = v[1];
                            up[i,j] = v[2];
                            # no down movement
                            southwest[i,j] = v[3];
                            southeast[i,j] = v[4]
                            northwest[i,j] = v[5];
                            northeast[i,j] = v[6];
                        end
    
                    # Now top edge
                    elseif 1<i<acrm && j == upn
                        if i%2 == 1 # Odd i index hexagons are different to even i index hexagons
                            v = rand(Multinomial(A[i,j],[1-(5/6)*pmig;pmig/6;pmig/6;pmig/6;pmig/6;pmig/6]));

                            remain[i,j] = v[1];
                            # no up movement
                            down[i,j] = v[2];
                            southwest[i,j] = v[3]
                            southeast[i,j] = v[4]
                            northwest[i,j] = v[5];
                            northeast[i,j] = v[6];

                        else
                            v = rand(Multinomial(A[i,j],[1-(3/6)*pmig;pmig/6;pmig/6;pmig/6]));
                            remain[i,j] = v[1];
                            # no up movement
                            down[i,j] = v[2]
                            southwest[i,j] = v[3];
                            southeast[i,j] = v[4]
                            # no northwest movement
                            # no northeast movement
                        end
    
                    elseif i == 1 && j == 1
                        #println("γ")
                        v = rand(Multinomial(A[i,j],[1-(3/6)*pmig;pmig/6;pmig/6;pmig/6]));
                        #println("δ")
                        remain[i,j] = v[1];
                        up[i,j] = v[2];
                        northeast[i,j] = v[3];
    
                    elseif i == acrm && j == 1
                        if i%2 == 1
                            v = rand(Multinomial(A[i,j],[1-(2/6)*pmig;pmig/6;pmig/6]));
                            remain[i,j] = v[1];
                            up[i,j] = v[2];
                            northwest[i,j] = v[3];
                            
                        else
                            v = rand(Multinomial(A[i,j],[1-(3/6)*pmig;pmig/6;pmig/6;pmig/6]));
                            remain[i,j] = v[1];
                            up[i,j] = v[2];
                            southwest[i,j] = v[3];
                            northwest[i,j] = v[4];
                        end
    
                    elseif i == 1 && j == upn
                        v = rand(Multinomial(A[i,j],[1-(3/6)*pmig;pmig/6;pmig/6;pmig/6]));
                        remain[i,j] = v[1];
                        down[i,j] = v[2];
                        southeast[i,j] = v[3];
                        northeast[i,j] = v[4];
    
                    elseif i == acrm && j == upn
                        if i%2 == 1
                            v = rand(Multinomial(A[i,j],[1-(3/6)*pmig;pmig/6;pmig/6;pmig/6]));
                            remain[i,j] = v[1];
                            down[i,j] = v[2];
                            southwest[i,j] = v[3];
                            northwest[i,j] = v[4];

                        else
                            v = rand(Multinomial(A[i,j],[1-(2/6)*pmig;pmig/6;pmig/6]));
                            remain[i,j] = v[1];
                            down[i,j] = v[2];
                            southwest[i,j] = v[3];
                        end
                    else
                        print("Error calculating emigrant numbers")
                        return 2
                    end
                
            end
        end
    
        #for i = 1:acrm
        #    if down[i,1] != 0
         #       print("Error calculating emigrant numbers")
         #       return 3
          #  elseif up[i,upn] != 0
         #       print("Error calculating emigrant numbers")
          #      return 4
        #    end
        #end
    
        #for j = 1:upn
        #    if left[1,j] != 0
        #        print("Error calculating emigrant numbers")
        #        return 5
        #    elseif right[acrm,j] != 0
        #        print("Error calculating emigrant numbers")
        #        return 6
        #    end
        #end
    
        # Now calculate population distribution after migration, B[i,j]

        uptemp = zeros(Int,acrm+2,upn+2);
        downtemp = zeros(Int,acrm+2,upn+2);
        southwesttemp = zeros(Int,acrm+2,upn+2);
        southeasttemp = zeros(Int,acrm+2,upn+2);
        northwesttemp = zeros(Int,acrm+2,upn+2);
        northeasttemp = zeros(Int,acrm+2,upn+2);

        uptemp[2:acrm+1,2:upn+1] = up;
        downtemp[2:acrm+1,2:upn+1] = down;
        southwesttemp[2:acrm+1,2:upn+1] = southwest;
        southeasttemp[2:acrm+1,2:upn+1] = southeast;
        northwesttemp[2:acrm+1,2:upn+1] = northwest;
        northeasttemp[2:acrm+1,2:upn+1] = northeast;

        B = zeros(Int,(acrm,upn));
    
        for i = 1:acrm
            for j = 1:upn
                if i%2 == 1
                    B[i,j] = remain[i,j] + uptemp[i+1,j] + downtemp[i+1,j+2] + southwesttemp[i+2,j+1] + southeasttemp[i,j+1] + northwesttemp[i+2,j] + northeasttemp[i,j];
                else
                    B[i,j] = remain[i,j] + uptemp[i+1,j] + downtemp[i+1,j+2] + southwesttemp[i+2,j+2] + southeasttemp[i,j+2] + northwesttemp[i+2,j+1] + northeasttemp[i,j+1];
                end
            end
        end
    
        return B
    
end

function sqmignewfrequency(Rm,α,xold,yold,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig)
    # Updates haplotype frequency in next generation for each (square) deme
    # Reproduction occurs in each deme, and then the progeny migrate

    (acrm,upn) = size(Nm);
    #println(1)

    xnew = zeros(acrm,upn,na);
    ynew = zeros(acrm,upn,na);
    Nfnew = zeros(Int128,(acrm,upn));
    Nmalenew = zeros(Int128,(acrm,upn));
    Nnew = zeros(Int128,(acrm,upn));
    #println(2)

    for iii = 1:acrm
        for jjj = 1:upn
            (xnew[iii,jjj,:],ynew[iii,jjj,:],Nfnew[iii,jjj],Nmalenew[iii,jjj],Nnew[iii,jjj]) = newfrequency(Rm,α,xold[iii,jjj,:],yold[iii,jjj,:],μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf[iii,jjj],Nm[iii,jjj],Gauss);
        end
    end
    #print(3)

    haplotype_number_x_new = zeros(Int128,(acrm,upn,na));
    haplotype_number_y_new = zeros(Int128,(acrm,upn,na));
    #print(4)

    for iii = 1:acrm
        for jjj = 1:upn
            for kkkk = 1:na
                # Check the following few lines!!!
                if isnan(xnew[iii,jjj,kkkk]*Nfnew[iii,jjj]) == 1
                    println(xold[iii,jjj,:])
                    println(yold[iii,jjj,:])
                    println(Nf[iii,jjj])
                    println(Nm[iii,jjj])
                    print(xnew[iii,jjj,:])
                    println(Nfnew[iii,jjj])
                    print(ynew[iii,jjj,:])
                    println(Nmalenew[iii,jjj])
                    #print("iii=")
                    #println(iii)
                    #print("jjj=")
                    #println(jjj)
                    #print("allele=")
                    #print(kkkk)
                    haplotype_number_x_new[iii,jjj,kkkk] = 0;
                    return (55,iii,jjj,kkkk)
                end
                if isnan(ynew[iii,jjj,kkkk]*Nmalenew[iii,jjj]) == 1
                    println(xold[iii,jjj,:])
                    println(yold[iii,jjj,:])
                    println(Nf[iii,jjj])
                    println(Nm[iii,jjj])
                    print(xnew[iii,jjj,:])
                    println(Nfnew[iii,jjj])
                    print(ynew[iii,jjj,:])
                    println(Nmalenew[iii,jjj])
                    println("allele=")
                    println(kkkk)
                    haplotype_number_y_new[iii,jjj,kkkk] = 0;
                    return (56,iii,jjj,kkkk)
                end
                haplotype_number_x_new[iii,jjj,kkkk] = round(Int128,xnew[iii,jjj,kkkk]*Nfnew[iii,jjj]);
                haplotype_number_y_new[iii,jjj,kkkk] = round(Int128,ynew[iii,jjj,kkkk]*Nmalenew[iii,jjj]);
            end
        end
    end
    #print(5)

    # Now migrate the haplotypes independently
    B_x = zeros(Int128,(acrm,upn,na));
    B_y = zeros(Int128,(acrm,upn,na));
    #println(6)
    for kkkk = 1:na
        A_x = haplotype_number_x_new[:,:,kkkk];
        A_y = haplotype_number_y_new[:,:,kkkk];
        #println(6.5)

        B_x[:,:,kkkk] = square_deme_migration(A_x,pmig,Gauss);
        B_y[:,:,kkkk] = square_deme_migration(A_y,pmig,Gauss);
    end
    #print(7)

    # Finally, calculate new population sizes and haplotype frequencies
    Nf_post_mig = zeros(Int128,(acrm,upn));
    Nmale_post_mig = zeros(Int128,(acrm,upn));
    N_post_mig = zeros(Int128,(acrm,upn))
    x_post_mig = zeros(acrm,upn,na);
    y_post_mig = zeros(acrm,upn,na);
    #print(8)

    for iii = 1:acrm
        for jjj = 1:upn
            Nf_post_mig[iii,jjj] = sum(B_x[iii,jjj,:]);
            Nmale_post_mig[iii,jjj] = sum(B_y[iii,jjj,:]);
            N_post_mig[iii,jjj] = Nf_post_mig[iii,jjj] + Nmale_post_mig[iii,jjj];
            if Nf_post_mig[iii,jjj] > 0
                x_post_mig[iii,jjj,:] = B_x[iii,jjj,:]./Nf_post_mig[iii,jjj];
            else
                x_post_mig[iii,jjj,:] = B_x[iii,jjj,:];
            end

            if Nmale_post_mig[iii,jjj] > 0
                y_post_mig[iii,jjj,:] = B_y[iii,jjj,:]./Nmale_post_mig[iii,jjj];
            else
                y_post_mig[iii,jjj,:] = B_y[iii,jjj,:];
            end
        end
    end
    #print(9)

    return (x_post_mig,y_post_mig,Nf_post_mig,Nmale_post_mig,N_post_mig,B_x,B_y)

end

function hexmignewfrequency(Rm,α,xold,yold,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig)
    # Updates haplotype frequency in next generation for each (square) deme
    # Reproduction occurs in each deme, and then the progeny migrate

    (acrm,upn) = size(Nm);
    #println(1)

    xnew = zeros(acrm,upn,na);
    ynew = zeros(acrm,upn,na);
    Nfnew = zeros(Int128,(acrm,upn));
    Nmalenew = zeros(Int128,(acrm,upn));
    Nnew = zeros(Int128,(acrm,upn));
    #println(2)

    for iii = 1:acrm
        for jjj = 1:upn
            (xnew[iii,jjj,:],ynew[iii,jjj,:],Nfnew[iii,jjj],Nmalenew[iii,jjj],Nnew[iii,jjj]) = newfrequency(Rm,α,xold[iii,jjj,:],yold[iii,jjj,:],μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf[iii,jjj],Nm[iii,jjj],Gauss);
        end
    end
    #print(3)

    haplotype_number_x_new = zeros(Int128,(acrm,upn,na));
    haplotype_number_y_new = zeros(Int128,(acrm,upn,na));
    #print(4)

    for iii = 1:acrm
        for jjj = 1:upn
            for kkkk = 1:na
                # Check the following few lines!!!
                if isnan(xnew[iii,jjj,kkkk]*Nfnew[iii,jjj]) == 1
                    println(xold[iii,jjj,:])
                    println(yold[iii,jjj,:])
                    println(Nf[iii,jjj])
                    println(Nm[iii,jjj])
                    print(xnew[iii,jjj,:])
                    println(Nfnew[iii,jjj])
                    print(ynew[iii,jjj,:])
                    println(Nmalenew[iii,jjj])
                    #print("iii=")
                    #println(iii)
                    #print("jjj=")
                    #println(jjj)
                    #print("allele=")
                    #print(kkkk)
                    haplotype_number_x_new[iii,jjj,kkkk] = 0;
                    return (55,iii,jjj,kkkk)
                end
                if isnan(ynew[iii,jjj,kkkk]*Nmalenew[iii,jjj]) == 1
                    println(xold[iii,jjj,:])
                    println(yold[iii,jjj,:])
                    println(Nf[iii,jjj])
                    println(Nm[iii,jjj])
                    print(xnew[iii,jjj,:])
                    println(Nfnew[iii,jjj])
                    print(ynew[iii,jjj,:])
                    println(Nmalenew[iii,jjj])
                    println("allele=")
                    println(kkkk)
                    haplotype_number_y_new[iii,jjj,kkkk] = 0;
                    return (56,iii,jjj,kkkk)
                end
                haplotype_number_x_new[iii,jjj,kkkk] = round(Int128,xnew[iii,jjj,kkkk]*Nfnew[iii,jjj]);
                haplotype_number_y_new[iii,jjj,kkkk] = round(Int128,ynew[iii,jjj,kkkk]*Nmalenew[iii,jjj]);
            end
        end
    end
    #print(5)

    # Now migrate the haplotypes independently
    B_x = zeros(Int128,(acrm,upn,na));
    B_y = zeros(Int128,(acrm,upn,na));
    #println(6)
    for kkkk = 1:na
        A_x = haplotype_number_x_new[:,:,kkkk];
        A_y = haplotype_number_y_new[:,:,kkkk];
        #println(6.5)

        B_x[:,:,kkkk] = hex_deme_migration(A_x,pmig);
        B_y[:,:,kkkk] = hex_deme_migration(A_y,pmig);
    end
    #print(7)

    # Finally, calculate new population sizes and haplotype frequencies
    Nf_post_mig = zeros(Int128,(acrm,upn));
    Nmale_post_mig = zeros(Int128,(acrm,upn));
    N_post_mig = zeros(Int128,(acrm,upn))
    x_post_mig = zeros(acrm,upn,na);
    y_post_mig = zeros(acrm,upn,na);
    #print(8)

    for iii = 1:acrm
        for jjj = 1:upn
            Nf_post_mig[iii,jjj] = sum(B_x[iii,jjj,:]);
            Nmale_post_mig[iii,jjj] = sum(B_y[iii,jjj,:]);
            N_post_mig[iii,jjj] = Nf_post_mig[iii,jjj] + Nmale_post_mig[iii,jjj];
            if Nf_post_mig[iii,jjj] > 0
                x_post_mig[iii,jjj,:] = B_x[iii,jjj,:]./Nf_post_mig[iii,jjj];
            else
                x_post_mig[iii,jjj,:] = B_x[iii,jjj,:];
            end

            if Nmale_post_mig[iii,jjj] > 0
                y_post_mig[iii,jjj,:] = B_y[iii,jjj,:]./Nmale_post_mig[iii,jjj];
            else
                y_post_mig[iii,jjj,:] = B_y[iii,jjj,:];
            end
        end
    end
    #print(9)

    return (x_post_mig,y_post_mig,Nf_post_mig,Nmale_post_mig,N_post_mig,B_x,B_y)

end

function fixedmig_TSR(m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,width,height,pmig,preexisting,T,test_threshold,Gauss,terminate_resistance)
    # m = number of target sites
    # s = homozygous fitness cost of drive
    # s*h = fitness cost of heterozygous W/D individual
    # s*hN = fitness cost of heterozygous W/N or R/D individual
    # σ_i = fitness cost of R allele at target site i, 1≦i≦m
    # Rm =
    # K =
    # xD,yD = initial proportion of D...D haplotypes in females/males respectively
    # ϵ = efficiency of drive cleavage per target site

    # ν = non-homologous end-joining (NHEJ) rate per generation per individual per
    # target site
    # μ = mutation rate per cleavage target site per generation per individual per
    # target site (e.g. if a target site contains 10bp and bp mutation rate is 1e-9,
    # μ=1e-8)
    # β = fraction of functional resistant NHEJ mutations
    # ξ = fraction of functional resistant single nucleotide mutations

    #0.05, 5.4e-8, 1e-4, 1,

    # preexisting = 1,0 to run a burn in period of 1/σ and set initial allele
    # frequencies or not, respectively
    # T = length of simulations in generations
    # test_threshold = threshold for the proportion of resistant haplotypes present
    # for resistance to be established
    # Gauss = 0 uses inbuilt multinomial generator; otherwise uses Gauss-Poisson
    # approximation
    # D = diffusivity of Brownian motion migration
    # pmig = probability opf migrating to an adjacent deme
    # area_approx = 1 approximates the deme size by calculating the size of the
    # equivalent circular deme
    
        #Set proportion of males
        ψ = 0.5;
    
        N = round(Int128,K);
    
        # To fit the area to square demes, we use the smallest number of demes in
        # each direction needed to cover the space
        acrm = width;
        upn = height;
        #println("checkpoint 1")
        A = zeros(Int128,(acrm,upn));
    
        # Now we distribute the individuals amongst the demes.
        # For now, assume uniform distribution. This can be changed later.
        K_perdeme = div(K,acrm*upn,RoundUp);
        if Gauss == 0 && log2(K_perdeme) > 63
            print("Error - if population > 2^63 then must use Gauss-Poisson Approximation")
            return -2
        end
        for i = 1:acrm
            for j = 1:upn
                A[i,j] = K_perdeme;
            end
        end
        #println("checkpoint 2")
    
        # Sort population per deme into males and females
        Nm = zeros(Int128,(acrm,upn));
        Nf = zeros(Int128,(acrm,upn));
    
        if Gauss == 0
            for i = 1:acrm
                for j = 1:upn
                    Nm[i,j] = rand(Binomial(A[i,j],ψ));
                    Nf[i,j] = A[i,j] - Nm[i,j];
                end
            end
        else
            for i = 1:acrm
                for j = 1:upn
                    vvvv = GaussPoissonHybrid_mnrnd(A[i,j],[ψ;1-ψ]);
                    Nm[i,j] = round(Int128(vvvv[1]));
                    Nf[i,j] = A[i,j] - Nm[i,j];
                end
            end
        end
    
        N = Nf + Nm;
        #println("checkpoint 3")
    
        # Density dependent parameter from Beverton-Holt model of population
        # dynamics
        # Note that we use the carrying capacity of *each deme*
        α = K_perdeme/(Rm - 1);
    
        # Dominance coefficients
        hD = h;
    
        if length(σ) == 1
            σ = σ*fill(1.0,(m,1));
        elseif length(σ) != m
            print("Error - σ should be a scalar of vector of length m")
            return -3
        end
        #println("checkpoint 4")
    
        # Set number of alleles (this would be difficult to change anyway in
        # practice)
        n = 4;
    
        # calculate number of (ordered) haplotypes
        na  = (n-1)^m + 1;
    
        #calculate fitness matrices
        Wm = fill(1.0,(na,na));
        Wf = femalefitnessmatrix(m,n,s,hD,hN,σ);
    
        #Calculate mutation probaility matrix. M[i,j] is the probability of de novo
        #mutation haplotype j -> haplotype i
        mutmatrix = mutationmatrix(m,n,μ,ξ);
    
        #Calculate drive conversion matrix. Suppose an individual has genotype
        #[j,drive]. Then K[i,j]+δ_{ij} is the probability of the conversion j->i via
        #drive and NHEJ mutations. So K[i,j] is the mean relative change in
        #frequency j->i.
        κ = drivematrix(m,n,ϵ,ν,β);
        #println("checkpoint 5")
    
        ## Initial frequencies
        x0 = zeros((acrm,upn,na));
        y0 = zeros((acrm,upn,na));
        for i = 1:acrm
            for j = 1:upn
                x0[i,j,1] = 1.0;
                y0[i,j,1] = 1.0;
            end
        end
        #println("checkpoint 6")
    
        ## If preexisting
        if preexisting == 1
            #Calculate burn in time
            if minimum(σ) <= 0
                print("Error - for preexisting = 1, resistance alleles must always incur fitness cost")
                return
            end
    
            Tburn = round(1/minimum(σ));
    
            # Now initialise as in original code
            xxx0 = zeros(na);
            yyy0 = zeros(na);
            xxx0[1] = 1.0;
            yyy0[1] = 1.0;
            #println("checkpoint 7")
    
            indm = mutmatrix.>0;
            for i = 1:na
                for j = 2:na
                    indm[i,j] = 0;
                end
            end
            indm[1,1] = 0;
    
            muin = mutmatrix[indm];
            #print("muin = ")
            #println(muin)
    
            ffff = diagm((Wf[indm].-1)./2);
            #print("ffff = ")
            #println(ffff)
            #println(ffff)
    
            indm2 = mutmatrix.>0;
            indm2[1,1] = 0;
            for i = 2:na
                for j = 2:na
                    if indm2[i,1] > 0 && indm2[j,1] > 0
                        indm2[i,j] = 1;
                    else
                        indm2[i,j] = 0;
                    end
                end
            end
    
            for i = 1:na
                indm2[i,1] = 0;
                indm2[1,i] = 0;
            end
    
            modmut1 = mutmatrix[indm2];
            lengthvec = length(modmut1);
            if lengthvec != length(ffff);
                println("Error with preexisting")
                println(lengthvec)
                println(length(ffff))
            end
            lengthvec = round(Int,sqrt(lengthvec));
            modmut2 = zeros((lengthvec,lengthvec));
            for i = 1:lengthvec
                for j = 1:lengthvec
                    modmut2[i,j] = modmut1[lengthvec*(j-1)+i];
                    if i == j
                        modmut2[i,j] += -1;
                    end
                end
            end
            #print("modmut2 = ")
            #println(modmut2)
    
            indmvec = BitArray(undef,na);
            for i = 1:na
                indmvec[i] = indm[i,1];
            end
    
            xxx0[indmvec] = -inv(modmut2+ffff)*muin;
            yyy0[indmvec] = -inv(modmut2+ffff)*muin;
            #println(x0)
            #println(y0)
    
            xxx0[na] = 0;
            yyy0[na] = 0;
            xxx0[1] = 1 - sum(xxx0[2:end]);
            yyy0[1] = 1 - sum(yyy0[2:end]);
            #println("checkpoint 8")
            print("xxx0 = ")
            println(xxx0)
    
            for i = 1:acrm
                for j = 1:upn
                    x0[i,j,:] = xxx0;
                    y0[i,j,:] = yyy0;
                end
            end
            #println("checkpoint 9, burn-in commencing")
    
            for tt = 1:Tburn
                println(tt)
                (x0,y0,Nf,Nm,N,~,~) = sqmignewfrequency(Rm,α,x0,y0,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig);
            end
            #println("checkpoint 10, burn-in complete")
        end
    
        ## Now run simulation, keeping track of frequencies and population size
        # First, add drive to initial frequencies
        #println("checkpoint 11")
        for i = 1:acrm
            for j = 1:upn
                x0[i,j,na] = xD;
                x0[i,j,1] = x0[i,j,1] - xD;
                y0[i,j,na] = yD;
                y0[i,j,1] = y0[i,j,1] - yD;
    
                if x0[i,j,1] < 0 || y0[i,j,1] < 0
                    println("Error - not enough wild alleles to add drive in deme")
                    println([i,j])
                    return 5
                end
            end
        end
        #println("checkpoint 12")
    
        xcurrent = x0;
        ycurrent = y0;
        #println("checkpoint 13")
    
        # x,y will keep track over time of fe/male haplotype frequencies,
        # respectively
        x = zeros(acrm,upn,na,T);
        y = zeros(acrm,upn,na,T);
        xnum = zeros(Int128, acrm,upn,na,T);
        ynum = zeros(Int128, acrm,upn,na,T);
        time = zeros(Int,T);
        x[:,:,:,1] = xcurrent;
        y[:,:,:,1] = ycurrent;
        #println("checkpoint 14")
    
        # popsize keeps track of female, male and total populations
        popsize = zeros(Int128,(3,acrm,upn,T));
        popsize[1,:,:,1] = Nf;
        popsize[2,:,:,1] = Nm;
        popsize[3,:,:,1] = Nf+Nm;
        #println("checkpoint 15")
    
        # We want to keep track of the resistant allele
        testallele = 2;
        resistance = 0;
        resistancepositions = resistancehaplotypes(m,n,[2]);
        #println("checkpoint 16")
    
        for t = 1:T-1
            #println(t)
            time[t+1] = t;
            #(xcurrent,ycurrent,Nf,Nm,N) = sqmignewfrequency(Rm,α,xcurrent,ycurrent,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig);
            newfrequencies = sqmignewfrequency(Rm,α,xcurrent,ycurrent,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig);
            if newfrequencies[1] == 55 || newfrequencies[1] == 56
                (~,iii,jjj,kkkk) = newfrequencies;
                println(x[iii,jjj,:,1:t+1])
                println(y[iii,jjj,:,1:t+1])
                println(popsize[:,iii,jjj,1:t+1])
                return (x[iii,jjj,:,1:t+1],y[iii,jjj,:,1:t+1],popsize[:,iii,jjj,1:t+1],resistance,τ)
            else
                (xcurrent,ycurrent,Nf,Nm,N,xnum_current,ynum_current) = newfrequencies;
            end
    
            x[:,:,:,t+1] = xcurrent;
            y[:,:,:,t+1] = ycurrent;
            xnum[:,:,:,t+1] = xnum_current;
            ynum[:,:,:,t+1] = ynum_current;
            popsize[1,:,:,t+1] = Nf;
            popsize[2,:,:,t+1] = Nm;
            popsize[3,:,:,t+1] = Nf+Nm;
    
            # Test for resistance
            # Discuss exactly how do do this with supervisor. For now, do the following
            resistfreqf = zeros(acrm,upn);
            resistfreqm = zeros(acrm,upn);
            for i = 1:acrm
                for j = 1:upn
                    resistfreqf[i,j] = dot(xcurrent[i,j,:],resistancepositions);
                    resistfreqm[i,j] = dot(ycurrent[i,j,:],resistancepositions);
                end
            end
            #scalarresistf = mean(resistfreqf);
            #scalarresistm = mean(resistfreqm);
            scalarresistf = 0;
            scalarresistm = 0;
            for i = 1:acrm
                for j = 1:upn
                    scalarresistf += resistfreqf[i,j]*Nf[i,j];
                    scalarresistm += resistfreqm[i,j]*Nm[i,j];
                end
            end
            scalarresistf = scalarresistf;
            scalarresistm = scalarresistm;
            if scalarresistf >= test_threshold  && scalarresistm >= test_threshold
                if resistance == 0
                    τ = t;
                end
                resistance = 1;
                if terminate_resistance == 1
                    return (x[:,:,:,1:t+1],y[:,:,:,1:t+1],popsize[:,:,:,1:t+1],resistance,τ)
                end
            end
    
            if maximum(Nf + Nm) == 0
                #println("Population zero at t = $(t+1)")
                break
            end
        end
        #println("checkpoint 17")
        if resistance == 0
            τ = NaN;
        end
    
        return (x,y,popsize,resistance,τ)
    
end

function hexmig_TSR(m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,width,height,pmig,preexisting,T,test_threshold,Gauss,terminate_resistance)
    # m = number of target sites
    # s = homozygous fitness cost of drive
    # s*h = fitness cost of heterozygous W/D individual
    # s*hN = fitness cost of heterozygous W/N or R/D individual
    # σ_i = fitness cost of R allele at target site i, 1≦i≦m
    # Rm =
    # K =
    # xD,yD = initial proportion of D...D haplotypes in females/males respectively
    # ϵ = efficiency of drive cleavage per target site
    # ν = non-homologous end-joining (NHEJ) rate per generation per individual per
    # target site
    # μ = mutation rate per cleavage target site per generation per individual per
    # target site (e.g. if a target site contains 10bp and bp mutation rate is 1e-9,
    # μ=1e-8)
    # β = fraction of functional resistant NHEJ mutations
    # ξ = fraction of functional resistant single nucleotide mutations
    # preexisting = 1,0 to run a burn in period of 1/σ and set initial allele
    # frequencies or not, respectively
    # T = length of simulations in generations
    # test_threshold = threshold for the proportion of resistant haplotypes present
    # for resistance to be established
    # Gauss = 0 uses inbuilt multinomial generator; otherwise uses Gauss-Poisson
    # approximation
    # D = diffusivity of Brownian motion migration
    # pmig = probability opf migrating to an adjacent deme
    # area_approx = 1 approximates the deme size by calculating the size of the
    # equivalent circular deme
    
        #Set proportion of males
        ψ = 0.5;
    
        N = round(Int128,K);
    
        # To fit the area to square demes, we use the smallest number of demes in
        # each direction needed to cover the space
        acrm = width;
        upn = height;
        #println("checkpoint 1")
        A = zeros(Int128,(acrm,upn));
    
        # Now we distribute the individuals amongst the demes.
        # For now, assume uniform distribution. This can be changed later.
        K_perdeme = div(K,acrm*upn,RoundUp);
        if Gauss == 0 && log2(K_perdeme) > 63
            print("Error - if population > 2^63 then must use Gauss-Poisson Approximation")
            return -2
        end
        for i = 1:acrm
            for j = 1:upn
                A[i,j] = K_perdeme;
            end
        end
        #println("checkpoint 2")
    
        # Sort population per deme into males and females
        Nm = zeros(Int128,(acrm,upn));
        Nf = zeros(Int128,(acrm,upn));
    
        if Gauss == 0
            for i = 1:acrm
                for j = 1:upn
                    Nm[i,j] = rand(Binomial(A[i,j],ψ));
                    Nf[i,j] = A[i,j] - Nm[i,j];
                end
            end
        else
            for i = 1:acrm
                for j = 1:upn
                    vvvv = GaussPoissonHybrid_mnrnd(A[i,j],[ψ;1-ψ]);
                    Nm[i,j] = round(Int128(vvvv[1]));
                    Nf[i,j] = A[i,j] - Nm[i,j];
                end
            end
        end
    
        N = Nf + Nm;
        #println("checkpoint 3")
    
        # Density dependent parameter from Beverton-Holt model of population
        # dynamics
        # Note that we use the carrying capacity of *each deme*
        α = K_perdeme/(Rm - 1);
    
        # Dominance coefficients
        hD = h;
    
        if length(σ) == 1
            σ = σ*fill(1.0,(m,1));
        elseif length(σ) != m
            print("Error - σ should be a scalar of vector of length m")
            return -3
        end
        #println("checkpoint 4")
    
        # Set number of alleles (this would be difficult to change anyway in
        # practice)
        n = 4;
    
        # calculate number of (ordered) haplotypes
        na  = (n-1)^m + 1;
    
        #calculate fitness matrices
        Wm = fill(1.0,(na,na));
        Wf = femalefitnessmatrix(m,n,s,hD,hN,σ);
    
        #Calculate mutation probaility matrix. M[i,j] is the probability of de novo
        #mutation haplotype j -> haplotype i
        mutmatrix = mutationmatrix(m,n,μ,ξ);
    
        #Calculate drive conversion matrix. Suppose an individual has genotype
        #[j,drive]. Then K[i,j]+δ_{ij} is the probability of the conversion j->i via
        #drive and NHEJ mutations. So K[i,j] is the mean relative change in
        #frequency j->i.
        κ = drivematrix(m,n,ϵ,ν,β);
        #println("checkpoint 5")
    
        ## Initial frequencies
        x0 = zeros((acrm,upn,na));
        y0 = zeros((acrm,upn,na));
        for i = 1:acrm
            for j = 1:upn
                x0[i,j,1] = 1.0;
                y0[i,j,1] = 1.0;
            end
        end
        #println("checkpoint 6")
    
        ## If preexisting
        if preexisting == 1
            #Calculate burn in time
            if minimum(σ) <= 0
                print("Error - for preexisting = 1, resistance alleles must always incur fitness cost")
                return
            end
    
            Tburn = round(1/minimum(σ));
    
            # Now initialise as in original code
            xxx0 = zeros(na);
            yyy0 = zeros(na);
            xxx0[1] = 1.0;
            yyy0[1] = 1.0;
            #println("checkpoint 7")
    
            indm = mutmatrix.>0;
            for i = 1:na
                for j = 2:na
                    indm[i,j] = 0;
                end
            end
            indm[1,1] = 0;
    
            muin = mutmatrix[indm];
            #print("muin = ")
            #println(muin)
    
            ffff = diagm((Wf[indm].-1)./2);
            #print("ffff = ")
            #println(ffff)
            #println(ffff)
    
            indm2 = mutmatrix.>0;
            indm2[1,1] = 0;
            for i = 2:na
                for j = 2:na
                    if indm2[i,1] > 0 && indm2[j,1] > 0
                        indm2[i,j] = 1;
                    else
                        indm2[i,j] = 0;
                    end
                end
            end
    
            for i = 1:na
                indm2[i,1] = 0;
                indm2[1,i] = 0;
            end
    
            modmut1 = mutmatrix[indm2];
            lengthvec = length(modmut1);
            if lengthvec != length(ffff);
                println("Error with preexisting")
                println(lengthvec)
                println(length(ffff))
            end
            lengthvec = round(Int,sqrt(lengthvec));
            modmut2 = zeros((lengthvec,lengthvec));
            for i = 1:lengthvec
                for j = 1:lengthvec
                    modmut2[i,j] = modmut1[lengthvec*(j-1)+i];
                    if i == j
                        modmut2[i,j] += -1;
                    end
                end
            end
            #print("modmut2 = ")
            #println(modmut2)
    
            indmvec = BitArray(undef,na);
            for i = 1:na
                indmvec[i] = indm[i,1];
            end
    
            xxx0[indmvec] = -inv(modmut2+ffff)*muin;
            yyy0[indmvec] = -inv(modmut2+ffff)*muin;
            #println(x0)
            #println(y0)
    
            xxx0[na] = 0;
            yyy0[na] = 0;
            xxx0[1] = 1 - sum(xxx0[2:end]);
            yyy0[1] = 1 - sum(yyy0[2:end]);
            #println("checkpoint 8")
            #print("xxx0 = ")
            #println(xxx0)
    
            for i = 1:acrm
                for j = 1:upn
                    x0[i,j,:] = xxx0;
                    y0[i,j,:] = yyy0;
                end
            end
            #println("checkpoint 9, burn-in commencing")
    
            for tt = 1:Tburn
                #println(tt)
                (x0,y0,Nf,Nm,N,~,~) = hexmignewfrequency(Rm,α,x0,y0,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig);
            end
            #println("checkpoint 10, burn-in complete")
        end
    
        ## Now run simulation, keeping track of frequencies and population size
        # First, add drive to initial frequencies
        #println("checkpoint 11")
        for i = 1:trunc(Int,acrm/2)
            for j = 1:trunc(Int,upn/2)
                x0[i,j,na] = xD;
                x0[i,j,1] = x0[i,j,1] - xD;
                y0[i,j,na] = yD;
                y0[i,j,1] = y0[i,j,1] - yD;
    
                if x0[i,j,1] < 0 || y0[i,j,1] < 0
                    println("Error - not enough wild alleles to add drive in deme")
                    println([i,j])
                    return 5
                end
            end
        end
        #println("checkpoint 12")
    
        xcurrent = x0;
        ycurrent = y0;
        #println("checkpoint 13")
    
        # x,y will keep track over time of fe/male haplotype frequencies,
        # respectively
        x = zeros(acrm,upn,na,T);
        y = zeros(acrm,upn,na,T);
        xnum = zeros(Int128, acrm,upn,na,T);
        ynum = zeros(Int128, acrm,upn,na,T);
        time = zeros(Int,T);
        x[:,:,:,1] = xcurrent;
        y[:,:,:,1] = ycurrent;
        #println("checkpoint 14")
    
        # popsize keeps track of female, male and total populations
        popsize = zeros(Int128,(3,acrm,upn,T));
        popsize[1,:,:,1] = Nf;
        popsize[2,:,:,1] = Nm;
        popsize[3,:,:,1] = Nf+Nm;
        #println("checkpoint 15")
    
        # We want to keep track of the resistant allele
        testallele = 2;
        resistance = 0;
        resistancepositions = resistancehaplotypes(m,n,[2]);
        #println("checkpoint 16")
    
        for t = 1:T-1
            #println(t)
            #time[t+1] = t;
            #(xcurrent,ycurrent,Nf,Nm,N) = sqmignewfrequency(Rm,α,xcurrent,ycurrent,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig);
            newfrequencies = hexmignewfrequency(Rm,α,xcurrent,ycurrent,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig);
            if newfrequencies[1] == 55 || newfrequencies[1] == 56
                (~,iii,jjj,kkkk) = newfrequencies;
                #println(x[iii,jjj,:,1:t+1])
                #println(y[iii,jjj,:,1:t+1])
                #println(popsize[:,iii,jjj,1:t+1])
                return (x[iii,jjj,:,1:t+1],y[iii,jjj,:,1:t+1],popsize[:,iii,jjj,1:t+1],resistance,τ)
            else
                (xcurrent,ycurrent,Nf,Nm,N,xnum_current,ynum_current) = newfrequencies;
            end
    
            x[:,:,:,t+1] = xcurrent;
            y[:,:,:,t+1] = ycurrent;
            xnum[:,:,:,t+1] = xnum_current;
            ynum[:,:,:,t+1] = ynum_current;
            popsize[1,:,:,t+1] = Nf;
            popsize[2,:,:,t+1] = Nm;
            popsize[3,:,:,t+1] = Nf+Nm;
    
            # Test for resistance
            # Discuss exactly how do do this with supervisor. For now, do the following
            resistfreqf = zeros(acrm,upn);
            resistfreqm = zeros(acrm,upn);
            for i = 1:acrm
                for j = 1:upn
                    resistfreqf[i,j] = dot(xcurrent[i,j,:],resistancepositions);
                    resistfreqm[i,j] = dot(ycurrent[i,j,:],resistancepositions);
                end
            end
            #scalarresistf = mean(resistfreqf);
            #scalarresistm = mean(resistfreqm);
            scalarresistf = 0;
            scalarresistm = 0;
            for i = 1:acrm
                for j = 1:upn
                    scalarresistf += resistfreqf[i,j]*Nf[i,j];
                    scalarresistm += resistfreqm[i,j]*Nm[i,j];
                end
            end
            scalarresistf = scalarresistf;
            scalarresistm = scalarresistm;
            if scalarresistf >= test_threshold  && scalarresistm >= test_threshold
                if resistance == 0
                    τ = t;
                end
                resistance = 1;
                if terminate_resistance == 1
                    return (x[:,:,:,1:t+1],y[:,:,:,1:t+1],popsize[:,:,:,1:t+1],resistance,τ)
                end
            end
    
            if maximum(Nf + Nm) == 0
                #println("Population zero at t = $(t+1)")
                break
            end
        end
        #println("checkpoint 17")
        if resistance == 0
            τ = NaN;
        end
    
        return (x,y,popsize,resistance,τ)
    
end

function square_deme_size(D,pmig)
    # Approximates the size of deme needed to *roughly* simulate a Brownian motion
    # with diffusivity D and probability of migration pmig
    

    r = 2*sqrt(-D * log(pmig));
    x = sqrt(pi) * r;
    
    # Output the deme size length, x
    return x
end

function hex_deme_size(D,pmig)
    # Approximates the size of deme needed to *roughly* simulate a Brownian motion
    # with diffusivity D and probability of migration pmig
    
    r = 2*sqrt(-D * log(pmig));
    x = sqrt(pi*2/(3* sqrt(3))) * r;


        return x
end

function sqmig_TSR(m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,width,height,D,pmig,preexisting,T,test_threshold,Gauss,terminate_resistance)

    # Diffusion constant-based equivalent of fixedmig_TSR. D is diffusion constant. width, height are actual width and height of region
    deme_length = square_deme_size(D,pmig);

    # To fit the area to square demes, we use the smallest number of demes in
    # each direction needed to cover the space
    acrm = div(width,deme_length,RoundUp);
    acrm = round(Int,acrm);
    upn = div(height,deme_length,RoundUp);
    upn = round(Int,upn);

    (x,y,popsize,resistance,τ)  =  fixedmig_TSR(m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,acrm,upn,pmig,preexisting,T,test_threshold,Gauss,terminate_resistance);

    return (x,y,popsize,resistance,τ)
end
    
function hexmig_diff_TSR(m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,width,height,D,pmig,preexisting,T,test_threshold,Gauss,terminate_resistance)

    # Diffusion constant-based equivalent of fixedmig_TSR. D is diffusion constant. width, height are actual width and height of region
    deme_length = hex_deme_size(D,pmig);

    # To fit the area to square demes, we use the smallest number of demes in
    # each direction needed to cover the space
    acrm = div(width,deme_length,RoundUp);
    acrm = round(Int,acrm);
    upn = div(height,deme_length,RoundUp);
    upn = round(Int,upn);

    (x,y,popsize,resistance,τ)  =  hexmig_TSR(m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,acrm,upn,pmig,preexisting,T,test_threshold,Gauss,terminate_resistance);

    return (x,y,popsize,resistance,τ)
end
    
