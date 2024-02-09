## file containing functions regarding the value function interaction

## all functions in here take a vector pz as given that FOR NOW DEPENDS ONLY ON THE CURRENT LEVEL OF Z

function guess_V(param::model,pzvec::Vector)
    @unpack_model param
    # guess an na x nz x nk matrix
end

function V_update(param::model,pz::Vector,Vnext::Matrix)
    @unpack_model param
    # update V using utility function and the transition probabilities
    # setup Vnew, apol:
    Vnew = zeros(na,nz,nk)
    apol = zeros(na,nz,nk)
    # find the continuation values:
    EV = cont_val(param,Vnext)
    Vchoices = zeros(na)
    # now loop over all states today, find maximum values:
    for (ia,iz,ik) in stateind # LOOP OVER STATES today
        # for state today, setup a na' x 1 vector of Value function:
        coh = w + (pz(iz) + γgrid(iz))*agrid[ia]
        # compute value at all possible a' choices:
        for iaprime in 1:na
            Vchoices[iaprime] = u(coh - agrid[iaprime]) + β*EV[iaprime,iz,ik]
        end
        # find the maximum value and the corresponding a':
        vn_now,ap_now = findmax(Vchoices)
        Vnew[ia,iz,ik] = vn_now
        apol[ia,iz,ik] = agrid[ap_now]   
    end
    return Vnew,apol
end

function cont_val(param::model,Vnext::Matrix)
    @unpack_model param
    # update V using utility function and the transition probabilities
    cval = zeros(na,nz,nk)
    ## loop over all possible states today:
    # Assets tomorrow, expectation over z today (goes to h or l) -- k today (goes to 0 or k+1)
    for (ia,iz,ik) in stateind 
        cval[ia,iz,ik] = sum( γprob[iz,izprime]*θ*Vnext[ia,izprime,max[nk,ik+1]] + (1-θ)*Vnext[ia,izprime,1] for izprime in 1:nz)
    end
    return cval
end

