
"""
Basic SIR simulation
    Used to generate data.
"""
function SIR_sim(N,beta,gamma) 
    Z = [1,0] 
    t = 0.0
    dt = 1.0

    a = zeros(2)
    Z2_out = Int64[]
    Z1_out = Int64[]
    count = 0

    while Z[2] < Z[1]
        a[1] = beta*(N-Z[1])*(Z[1]-Z[2])/(N-1)
        a[2] = gamma*(Z[1]-Z[2])

        a0 = a[1]+a[2]

        t = t - log(rand())/a0

        while t > dt*count
            push!(Z1_out,Z[1])
            push!(Z2_out,Z[2])
            count = count + 1
        end

        if rand()*a0 < a[1]
            Z[1] = Z[1] + 1
        else
            Z[2] = Z[2] + 1
        end
    end

    push!(Z1_out,Z[1])
    push!(Z2_out,Z[2])

    return (Z1_out,Z2_out)
end


"""
Forward substitution with no Q matrix.
    This version iterates over the whole state space. 
    this works by calculting on the fly the indices of the elements that are updated
    rather than these being stored in a sparse matrix. It seems that in practice, this 
    recalulation is just as fast as looking up the elements from the sparse CSC structure. 
    This only works becasue Q has a lower triangular structure. 
"""
function fwd_noQ(phi,beta::Float64,gamma::Float64,N)

    for i=1:100
        count = 2;
        for z1 = 1:N-1   # dont do the first element as absorbing. 
            S = N-z1
            for z2 = 0:z1-1   # dont do top of column as absorbing

                a1 = beta*S*(z1-z2)
                a2 = gamma*(z1-z2)

                phi[count] = phi[count] / (1.0+ (a1+a2))
                phi[count+1] = phi[count+1] + phi[count]*a2
                phi[count+z1+1] = phi[count+z1+1] + phi[count]*a1
                count = count + 1
            end
            count = count + 1
        end
        # the last column which is a bit differnt. 
        z1 = N
        S = N-z1
        for z2 = 0:z1-1

            a2 = gamma*(z1-z2)

            phi[count] = phi[count] / (1+ a2 )
            phi[count+1] = phi[count+1] + phi[count]*a2
            count = count + 1
        end
    end 
end

"""
Forward substitution with no Q matrix.
    This version only iterates between given values of z1, which makes it a lot more
    efficient if we have exact observations of the state. 
"""
function fwd_noQ_tight(phi,beta::Float64,gamma::Float64,N,n1,n2)

    for i=1:100

        count = n1*(n1+1)÷2 +1
        
        for z1 = n1:n2-1
            S = N-z1
            for z2 = 0:z1-1   # dont do top of column as absorbing

                a1 = beta*S*(z1-z2)
                a2 = gamma*(z1-z2)

                phi[count] = phi[count] / (1+ (a1+a2))
                phi[count+1] = phi[count+1] + phi[count]*a2
                phi[count+z1+1] = phi[count+z1+1] + phi[count]*a1
                count = count + 1
            end
            count = count + 1
        end

        if n2 < N
            z1 = n2
            S = N-z1
            for z2 = 0:z1-1   # dont do top of column as absorbing

                a1 = beta*S*(z1-z2)
                a2 = gamma*(z1-z2)

                phi[count] = phi[count] / (1+ (a1+a2))
                phi[count+1] = phi[count+1] + phi[count]*a2
                phi[count+z1+1] = phi[count+z1+1] + phi[count]*a1
                count = count + 1
            end

        else
            # last col is a bit different. 
            z1 = N
            for z2 = 0:z1-1

                a2 = gamma*(z1-z2)

                phi[count] = phi[count] / (1+ a2)
                phi[count+1] = phi[count+1] + phi[count]*a2
                count = count + 1
            end
        end
    end 
end



"""
Forward substitution for final size calc.
    For doing the last part of the likelihood.
    phi will no longer be a prob vector after this.
"""
function fwd_fs_noQ(phi,beta::Float64,gamma::Float64,N)

    count = 2;
    for z1 = 1:N-1   # don't do the first element as absorbing. 
        S = N-z1
        for z2 = 0:z1-1   # dont do top of column as absorbing

            a1 = beta*S*(z1-z2)
            a2 = gamma*(z1-z2)
            p2 = a2/(a1+a2)

            phi[count+1] = phi[count+1] + phi[count]*p2
            phi[count+z1+1] = phi[count+z1+1] + phi[count]*(1-p2)
            count = count + 1
        end
        count = count + 1
    end

    # last colum is differnt as no infection events. 
    z1 = N
    S = N-z1
    for z2 = 0:z1-1

        phi[count+1] = phi[count+1] + phi[count]
        count = count + 1
    end

end



""""
Calculate the log-likelihood for the SIR model with infection events observed.
"""
function likelihood(beta::Float64, gamma::Float64, Z1_obs::Vector{Int64}, N::Int64)
  
    bigK = (N+1)*(N+2)÷2   
    steps = 100
    t_step = 1.0/ steps
    likes = zeros(length(Z1_obs) + 1)

    # initial condition is a single infected
    phi = zeros(bigK, 1)
    phi[2] = 1.0
  
    for i = 1:length(Z1_obs)

        fwd_noQ(phi,beta/(N-1)*t_step, gamma*t_step,N)

        # condition the prob vector
        c = (Z1_obs[i]+1)*Z1_obs[i]÷2
        for j=1:c
            phi[j] = 0.0
        end
        for j=c+Z1_obs[i]+2:bigK
            phi[j] = 0.0
        end

        # renormalise
        likes[i] = sum(phi)
        phi = phi / likes[i]

    end

    fwd_fs_noQ(phi,beta/(N-1),gamma,N)
    likes[end] = phi[(Z1_obs[end]+1)*(Z1_obs[end]+2)÷2]

    return sum(log.(likes))

  end

""""
Calculate likelihood.
  Efficient version.
  Calcualte the log-likelihood of a given sequence of Z1 observations.
"""
function likelihood_tight(beta::Float64, gamma::Float64, Z1_obs::Vector{Int64}, N::Int64)
     
    bigK = (N+1)*(N+2)÷2   
    steps = 100
    t_step = 1.0/ steps
    likes = zeros(length(Z1_obs)+1)

    # initial condition is a single infected
    phi = zeros(bigK, 1)
    phi[2] = 1.0

    # first step is different due to initial condition.
    fwd_noQ_tight(phi,beta/(N-1)*t_step, gamma*t_step,N,1,Z1_obs[1])

    c = (Z1_obs[1]+1)*Z1_obs[1]÷2
    for j=c+Z1_obs[1]+2:bigK
        phi[j] = 0.0
    end
    likes[1] = sum(phi)
    phi = phi / likes[1]
    
    for i = 2:length(Z1_obs)

        fwd_noQ_tight(phi,beta/(N-1)*t_step, gamma*t_step,N,Z1_obs[i-1],Z1_obs[i])

        # condition the prob vector
        c = (Z1_obs[i]+1)*Z1_obs[i]÷2
        for j=1:c
            phi[j] = 0.0
        end
        for j=c+Z1_obs[i]+2:bigK
            phi[j] = 0.0
        end

        # renormalise
        likes[i] = sum(phi)
        phi = phi / likes[i]

    end
  
      fwd_fs_noQ(phi,beta/(N-1),gamma,N)
      likes[end] = phi[(Z1_obs[end]+1)*(Z1_obs[end]+2)÷2]
  
      return sum(log.(likes))
  
    end