
"""
simulate forward 1 day until part+1 matches have been seen
"""
function forward_day_alive(N,beta,gamma,X,Z1_cum,y)

    part = size(X,1); #number of rows
    a = zeros(2);
    Z = zeros(Int,2)
    X1 = zeros(part); # holds the particles that match.
    K = 10^5;

    jj = 1;
    n = 0;    #% counts number of times we need to simulate until we get part+1 matches
    while jj <= part + 1
    
        # sample a particle
        ii = rand(1:part)
        n = n + 1;

        Z[1] = Z1_cum
        Z[2] = X[ii]; 
        t = 0.0;
        
        while Z[1] <= y && Z[1]-Z[2]> 0 

            I = Z[1]-Z[2];
            
            a[1] =  beta*(N-Z[1])*I/(N-1);
            a[2] =  gamma*I;
            a0 = a[1]+a[2];
  
            t = t - log(rand())/a0;
            
            if t > 1
                break
            end
                        
            if a0*rand() <= a[1]
                Z[1] = Z[1] + 1;
            else
                Z[2] = Z[2] + 1;
            end

        end
        
        if Z[1] == y && Z[1]-Z[2] > 0
            # particle matches observation
            # only add the first part matches and not the last.
            if jj <= part
                X1[jj] = Z[2];
            end
            jj = jj+1;
        end
              
        if n > K
            # stop it doing too many iterations and getting stuck
            return (0,X1);
        end

    end

    like = part/(n-1);
    
    return (like,X1)
end

"""
Estimate the prob of no further infections after last observation.
"""
function end_alive(N,beta,gamma,X,Z1_cum)


    part = size(X,1);
    a = zeros(2);
    Z = zeros(Int,2)
    K = 10^5;

    jj = 1;
    n = 0;  
    while jj <= part + 1

        ii = rand(1:part)
        n = n + 1;

        Z[1] = Z1_cum
        Z[2] = X[ii]; 
        t = 0.0;
        
        # no need to track time in this version.
        while Z[1] == Z1_cum && Z[1]-Z[2] > 0 

            I = Z[1]-Z[2];
            
            a[1] =  beta*(N-Z[1])*I/(N-1);
            a[2] =  gamma*I;
            a0 = a[1]+a[2];
                        
            if a0*rand() <= a[1]
                Z[1] = Z[1] + 1;
            else
                Z[2] = Z[2] + 1;
            end

        end
        
        if Z[1] == Z1_cum && Z[1]-Z[2] == 0
            jj = jj+1;
        end
              
        if n > K
            # stop it doing too many iterations and getting stuck
            return 0
        end

    end

    like = part/(n-1);
    
    return like
end

"""
Alive particle filter for basic SIR model.
"""
function SIR_likelihood_alive(N,beta,gamma,y,part)

    l = length(y)  # length of timeseries.
    LL = zeros(l)
    
    # set initial state, assume single infection event 
    Z2 = zeros(Int,part)
    # assume y is cumulative
    Z1 = 1

    for jj=1:l

        (like,Z2) = forward_day_alive(N,beta,gamma,Z2,Z1,y[jj])
        Z1 = y[jj]

        if like > 0
            LL[jj] = log(like)
        else
            return -Inf
        end   
    end
    
    # prob of no more events after the last (cumulative does not increase.)
    # don't need the particles after this step.
    LL[end] = log(end_alive(N,beta,gamma,Z2,y[end]))
    
    return sum(LL)
    
end

