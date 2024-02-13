# file containing all aggregate functions regarding asset demand, finding the price etc.

# function to guess the joint distribution:
function jdguess(param::model)
    @unpack_model param
    jd = zeros(na,nz,nk) # assets, last state observed, periods since info update
    # guess same weight for all states:
    for (ia,iz,ik) in stateind
        jd[ia,iz,ik] = 1/(na*nz*nk)
    end
    return jd
end

# function to move the jd one period forward:
function jdtrans(param::model,jd,apol,zid)
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
