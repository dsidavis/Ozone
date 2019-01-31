deltaX = 
function(UOS, n_t = length(UOS), fev_base, #input
               K # parameters
               )
    # Takes vectors of UOS and calculates FEV
{
#    r = (1 - exp(-K))
    ceoss = UOS / K

    .Call("R_deltaXLoop", ceoss, fev_base, K)
}    
