# Imprtance sampling particle filter for the SIR model
# Andrew Black - 25/06/2020

using StatsBase
using Distributions

"""
Forwards simulation using importance sampling
"""
function forward_day_is!(S,beta,gamma,X,Z1_cum,omega,NR)

    part = size(X,1); #number of rows
    a = zeros(2);
    Z = zeros(Int,2)
    
    for kk = 1:part

        Z[1] = Z1_cum
        Z[2] = X[kk]
        N = S
        
        # generate the infection times
        tR = sort!(rand(NR));

        r = 1;
        L_imp = 0 #-gammaln(NR+1);
        t = 0;
        
        while r <= NR

            I = Z[1]-Z[2];
            
            a[1] =  beta*(N-Z[1])*I;
            a[2] =  gamma*I
            a0 = a[1]+a[2];

            if I > 1
                b = a[2];
            else
                b = 0.0;
            end

            t_dash = -log(rand())/b;

               if t_dash + t < tR[r]
                   # recovery event

                    Z[2] = Z[2] +1;
                    t = t + t_dash;

                    L_imp = L_imp + log(a[2])-a0*(t_dash) - (log(b) - b*t_dash);
               else

                   L_imp = L_imp + log(a[1])-a0*(tR[r]-t) +b*(tR[r]-t);

                   Z[1] = Z[1] +1;
                   t = tR[r];     
                   r = r+1;

               end

        end

        # after the end of all the infection events then er need to put in the
        # last recoveries.
        # in this part we can't let the process fade out.

        while t < 1
             
            I = Z[1]-Z[2];
            a[1] =  beta*(N-Z[1])*I;
            a[2] =  gamma*I;
            a0 = a[1]+a[2];
 
            if I > 1
                b = a[2];
            else
                b = 0;
            end
            
            
            t_dash = -log(rand())/b;
            
            if t + t_dash > 1
                # gone over the next day
         
                L_imp = L_imp  -a0*(1-t) + b*(1-t);
                break;
                
            else
                # still within the day
                Z[2] = Z[2] +1;
                t = t + t_dash;

                L_imp = L_imp + log(a[2]) -a0*t_dash  - (log(b) - b*t_dash);
            end

        end

        omega[kk] = L_imp
        X[kk] = Z[2];      
    end
end


"""
Forwards simulation using importance sampling
"""
function forward_lastday_is!(S,beta,gamma,X,Z1_cum,omega,NR)

    part = size(X,1); #number of rows
    a = zeros(2);
    Z = zeros(Int,2)
    
    for kk = 1:part

        Z[1] = Z1_cum
        Z[2] = X[kk]
        N = S
        
        # generate the infection times
        tR = sort!(rand(NR));

        r = 1;
        L_imp = 0 #-gammaln(NR+1);
        t = 0;
        
        while r <= NR

            I = Z[1]-Z[2];
            
            a[1] =  beta*(N-Z[1])*I;
            a[2] =  gamma*I
            a0 = a[1]+a[2];

            if I > 1
                b = a[2];
            else
                b = 0.0;
            end

            t_dash = -log(rand())/b;

               if t_dash + t < tR[r]
                   # recovery event

                    Z[2] = Z[2] +1;
                    t = t + t_dash;

                    L_imp = L_imp + log(a[2])-a0*(t_dash) - (log(b) - b*t_dash);
               else

                   L_imp = L_imp + log(a[1])-a0*(tR[r]-t) +b*(tR[r]-t);

                   Z[1] = Z[1] +1;
                   t = tR[r];     
                   r = r+1;

               end

        end

        # as this is the last day,we can allow fadeout
        while t < 1
             
            I = Z[1]-Z[2];
            a[1] =  beta*(N-Z[1])*I;
            a[2] =  gamma*I;
            a0 = a[1]+a[2];
 
            b = a[2]

            t_dash = -log(rand())/b;
            
            if t + t_dash > 1
                # gone over the next day
         
                L_imp = L_imp  -a0*(1-t) + b*(1-t);
                break;
                
            else
                # still within the day
                Z[2] = Z[2] +1;
                t = t + t_dash;

                L_imp = L_imp + log(a[2]) -a0*t_dash  - (log(b) - b*t_dash);
            end

        end

        omega[kk] = L_imp
        X[kk] = Z[2];      
    end
end

"""
Return the probability of the epidemic fading out from the current state,
given no more infections occur. 
"""
function forward_end!(N,beta,gamma,X,Z1_final,omega)
    
    part = size(X,1)
    p = gamma/(beta*(N-Z1_final) + gamma)

    # no more infection events so Z1 does not change
    # which simplifies the calculation a lot. 

    for kk = 1:part
        omega[kk] =  (Z1_final-X[kk])*log(p)
    end
    
end
    
"""
Log-sum-exp trick.
weights are exponentiated in place ready for resampling.
"""
function log_sum_exp!(omega)

    os,_ = findmax(omega)

    for i=1:length(omega)
        omega[i] = exp(omega[i] - os)
    end    
    return log(sum(omega)) + os
end

"""
systematic resample. 
Done in place for the win. 
    pass sum_w from previous step to save recalculating this. 
    w - array of weights.
    X - particle array.
"""
function systematic_resample!(w,X)

    K = length(w)
    inc = sum(w) / K
    ac = inc * rand()

    i = 1; j=1

    upper_sum = w[1]

    while (i <= K)

        if (ac <= upper_sum)

            X[i] = X[j]
            ac = ac + inc
            i = i +1
        else
            j = j +1
            upper_sum += w[j]
        
        end
    end
end


"""
Particle filter for basic SIR model using importance sampling.
    In this version, y is a vector of daily counts, not a cumulative 
    number. Also it doesn't include the initial condition (might change
    at some point for flexibility).
"""
function SIR_likelihood_is(N,beta_un,gamma,y,part)
    
        l = length(y); 
     
        Z2 = zeros(Int,part)
        omega = zeros(part)

        # normalise here 
        beta = beta_un/(N-1);
        
        LL = zeros(l+1);
        
        #set initial state, single infected, no recovered.
        Z1 = 1
        
        for j=1:l-1
    
            forward_day_is!(N,beta,gamma,Z2,Z1,omega,y[j])
           
            # calcualte the log marginal likelihood
            # I only calculate the sum, not the mean, as this just translates 
            # into a constant factor in the likelihood.
            LL[j] = log_sum_exp!(omega)
            Z1 = Z1 + y[j]
    
            systematic_resample!(omega,Z2)
                
        end
        
        # on the last day the disease can fade-out after the last infection event.
        forward_lastday_is!(N,beta,gamma,Z2,Z1,omega,y[end])
        LL[l] = log_sum_exp!(omega)
        Z1 = Z1 + y[end]
        systematic_resample!(omega,Z2)

        # calculate the loglike of no futher infection events. 
        forward_end!(N,beta,gamma,Z2,Z1,omega)
        LL[l+1] = log_sum_exp!(omega)

        return sum(LL)
    
    end
    