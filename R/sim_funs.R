# Simulation of data with normal error for testing

sim_data = function(O3, Ve, CD, t, #inputs
                    Dos, K, A, n, # parameters
                    sigma #error )

{
    uos = UOS(O3, Ve, CD, t, Dos)
    x = FEV(uos, K, A)

    rnorm(n, x, sigma)       
}


