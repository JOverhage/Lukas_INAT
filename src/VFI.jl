## file containing functions regarding the value function interaction

## all functions in here take a vector pz as given that FOR NOW DEPENDS ONLY ON THE CURRENT LEVEL OF Z

function VFI(param::model,pzg)
    @unpack_model param
    # guess VF:
    Vg = guess_V(param,pzg)
    # setup matrices:
    Vnew = zeros(na,nz,nz,nk)
    apol = zeros(na,nz,nz,nk)
    iter = 0
    diff = 9999
    # iterate until convergence:
    for it in 1:maxit
        iter = iter + 1
        # update V using utility function and the transition probabilities
        Vnew,apol = V_Update(param,pzg,Vg)
        # compute diff:
        diff = maximum(abs.(Vnew - Vg))
        # print iteration and diff:
        println("Iteration: ",it," Diff: ",diff)
        # check convergence:
        if diff < tol
            # println("Value Function Converged, Iteration: ",it)
            break
        end
        # update guess:
        Vg = Vnew
    end
    # print iterations, diff:
    println(" -- VFI done, Iteration: ",iter, " Diff: ", diff)
    return Vnew,apol
end


function guess_V(param::model,pzvec::Vector)
    @unpack_model param
    # guess: Consume 10% of Assets, ignore any inattention and future values:
    Vg = zeros(na,nz,nz,nk)
    for (ia,iz,izl,ik) in stateind
        Vg[ia,iz,izl,ik] = u(param,w + (pzvec[iz] + γgrid[iz])*agrid[ia])
    end
    return Vg
end

function V_Update(param::model,pz::Vector,Vnext)
    @unpack_model param
    # update V using utility function and the transition probabilities
    # setup Vnew, apol:
    Vnew = zeros(na,nz,nz,nk)
    apol = zeros(na,nz,nz,nk)
    # find the continuation values:
    EV = cont_val(param,Vnext)
    Vchoices = zeros(na)
    # now loop over all states today, find maximum values:
    for (ia,iz,izl,ik) in stateind # LOOP OVER STATES today
        # for state today, setup a na' x 1 vector of Value function:
        coh = w + (pz[iz] + γgrid[iz])*agrid[ia]
        # compute value at all possible a' choices:
        for iaprime in 1:na
            Vchoices[iaprime] = u(param,coh - agrid[iaprime]) + β*EV[iaprime,iz,izl,ik]
        end
        # find the maximum value and the corresponding a':
        vn_now,ap_now = findmax(Vchoices)
        Vnew[ia,iz,iz,ik] = vn_now
        apol[ia,iz,iz,ik] = agrid[ap_now]   
    end
    return Vnew,apol
end

function cont_val(param::model,Vnext)
    @unpack_model param
    # update V using utility function and the transition probabilities
    cval = zeros(na,nz,nz,nk)
    ## loop over all possible states today:
    # Assets tomorrow, expectation over z today (goes to h or l) -- k today (goes to 0 or k+1)
    # transition matrix: take k period forward expectation (for kgrid = 0 just equal to true transition matrix), of going from izl (last observed state) to izprime (next state)
    for (ia,iz,izl,ik) in stateind 
        cval[ia,iz,izl,ik] = sum( γprob_kfwd[ik][izl,izprime]*θ*Vnext[ia,izprime,izl,min(nk,ik+1)] + γprob_kfwd[ik][izl,izprime]*(1-θ)*Vnext[ia,izprime,izprime,1] for izprime in 1:nz)
    end
    # transition matrix used is not 
    return cval
end

