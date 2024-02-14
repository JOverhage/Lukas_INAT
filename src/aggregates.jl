# file containing all aggregate functions regarding asset demand, finding the price etc.

# function to guess the joint distribution:
function jdguess(param::model)
    @unpack_model param
    jd = zeros(na,nz,nk) # assets, last state observed, periods since info update
    # guess same weight for all states:
    jd .= 1/(na*nz*nk)
    return jd
end

# function to move the jd one period forward:
function jdtrans(param::model,jd,apol,zid::Int64)
    @unpack_model param
    apnow = apol[:,:,zid,:] # evalate apolicy at the correct level of z today
    jdnext = copy.(jd) # copy the joint distribution
    for (ia,iz,ik) in stateind
        # Find the index of the point on agrid closest to current policy value:
        apind = findmin(abs.(agrid .- apnow[ia,iz,ik]))[2]
        # move the mass of non updaters to that asset grid point, at k' = k+1, last seen z unchanged:
        jdnext[apind,iz,ik+1] = θ*(jdnext[apind,iz,ik+1] + jd[ia,iz,ik])
        # move the mass of updaters to that asset grid point, at k'=0, last seen z = zid:
        jdnext[apind,zid,1] = (1-θ)*(jdnext[apind,iz,1] + jd[ia,iz,ik])
    end
    return jdnext
end

# function to compute total asset demand:
function asset_demand(param::model,jd,apol,znow::Int64)
    @unpack_model param
    # get correct policy:
    apnow = apol[:,:,znow,:] # evalate apolicy at the correct level of z today
    # create a mesh of asset demands given jd:
    admesh = apnow .* jd
    return sum(admesh)
end

#  price finder: takes guess as input, computes demand and bisects to find p
function p_finder(param::model,pgup,pgdown,zseries)
    @unpack_model param
    iter = 0
    # non-burn-in shocks:
    zT = zseries[Tburn+1:end]
    ## initial: Check bisection boundaries
    # simulate, and compute asset demand for pgup:
    jdup,apup,ADup = simulate(param,pgup,zseries) 
    # simulate, and compute asset demand for pgdown:
    jddown,apdown,ADdown = simulate(param,pgdown,zseries)
    # Intuition as usual High price => Low demand, low price => high demand, check if Dup < 1 < Ddown, if not print error:
    if (mean(ADup) < 1) || (mean(ADdown) > 1)
        println("Initial Bisection Boundaries not set correctly, ADup: ",mean(ADup)," ADdown: ",mean(ADdown))
        return pgup
    end
    # initialize loop midpoints:
    pmid = 0
    ADmid = 0
    ## Main loop:
    for it = 1:maxit_bis
        iter = iter + 1
        # compute midpoint:
        pmid = (pgup + pgdown)./2
        # simulate, and compute asset demand for pmid:
        jdmid,apmid,ADmid = simulate(param,pmid,zseries)
        # check if we are close enough:
        if sum(abs(ADmid[zT .== 1] - 1) + abs(ADmid[zT .== 2] - 1)) < tol_bis
            return pmid
        end
        # if not, update the boundaries for bad state:
        if ADmid[zT .== 1] > 1
            pgup[1] = pmid[1]
        else
            pgdown[1] = pmid[1]
        end
        # if not, update the boundaries for good state:
        if ADmid[zT .== 2] > 1
            pgup[2] = pmid[2]
        else
            pgdown[2] = pmid[2]
        end
    end
    printl(" -- Bisection ended at iteration: ", iter, " Excess demand g: ", mean(ADmid[zT .== 1] - 1), " Excess demand b: ", mean(ADmid[zT .== 2] - 1))
    return pmid
end

# simulation fct:
function simulate(param::model,pvec,zseries)
    @unpack_model param
    # simulate the model:
    # initialize the simulation:
    jdhist = [zeros(na,nz,nk) for i in 1:TT+Tburn]
    # initialize the joint distribution:
    jdhist[1] = jdguess(param)
    # for the given price, solve the VFI:
    V,ap = VFI(param,pvec)
    # loop over time:
    for t in 2:TT+Tburn
        # move the distribution forward, using the correct policy:
        jdhist[t] = jdtrans(param,jdhist[t-1],ap,zseries[t])
    end
    # compute asset demand:
    AD = [asset_demand(param,jdhist[t],ap,zseries[t]) for t in 1:TT+Tburn]
    # now return all time series beyond the burn in:
    return jdhist[Tburn+1:end],ap,AD[Tburn+1:end]
end

# drawing a series of z values:
function drawz(param::model)
    @unpack_model param
    # zinit = 1
    zdraws = [0 for i in 1:TT+Tburn]
    zdraws[1] = 1
    for t in 1:TT+Tburn-1
        # draw number between zero and one:
        draw = rand()
        # if draw is less than the transition probability, go to state 1, else go to state 2:
        if draw < γprob[zdraws[t]]
            zdraws[t+1] = 1
        else
            zdraws[t+1] = 2
        end
    end
    return zdraws
end



