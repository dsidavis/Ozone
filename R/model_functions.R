# Structural functions for each model

################################################################################
# From Schelegle et al.

UOS = function(O3, Ve, t, # inputs
               DR = O3 * Ve * 1.96,
               CD = cuml_integral(DR, t),
               Dos) #parameter
    # Vectorized means for calc UOS over a number of t timepoints
{
    DR / (1 + exp(-20 * (t - (Dos * t / CD))))
}

deltaX = function(UOS, n_t = length(UOS), fev_base, #input
               K # parameters
               )
    # Takes vectors of UOS and calculates FEV
{
    r = (1 - exp(-K))
    
    x = numeric(n_t)

    # Not right way to do this - every person starts with a baseline
    x[1] = fev_base
    

    for(i in 2:n_t)
        x[i] = x[i - 1] + (((UOS[i - 1] / K) - x[i - 1]) * r)

    
    cumsum(x)
}

FEV1 = function(dX, A){
    dX * A
}

experimentFEV1 = function(O3, Ve, t_stop, Dos, K, A)
    # O3, Ve are vectors with the measurement for the time interval
    # t_stop are the stop points (in min) for each associated interval, with t_stop[1] == 0
    # 

{
    dFEV1 = numeric(length(t_stop))
    dFEV1[1] = 0 # start at 0 delta FEV1
    
    for(i in seq(length(t_stop))) {
        t = t_stop[i]:t_stop[i+1]
        tmp = FEV(UOS(O3[i], Ve[i], t, Dos = dos), fev_base = dFEV[i-1])
        dFEV1[i] = FEV()
    }
    
                

}


cuml_integral = function(x, t)
    # Uses the trapzoid - integral approx. over a vector of x and t vals
{
    cumsum((c(0,x[-length(x)]) + x) / 2 * t)
}

################################################################################
# From McDonnell et al, 2013

FEV2 = function(Ui, M)
    # Non-linear form of the model
    # The E term is error, modeled elsewhere
{
    exp(Ui) * M
}


M = function(b1, b2, b3, b4, X, age, bmi)
    # Calculates the M values with covars
{
    (b1 + b2 * (age - 23.8) + b8 * (bmi - 23.1)) / (1 + b4 * exp(-b3 * X)) -
        (b1 + b2 * (age -23.8) + b8 * (bmi - 23.1)) / (1 + b4)
}

              
