################################################################################
# From McDonnell et al, 2013
if(FALSE){
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

              
}


XB5 = XB5 * ( exp( -B5 * TD ) ) +
    ( Cm * Vm * (B5^-1) ) * ( 1 - exp(-B5 * TD) ) +
    ( Cm * Vs * (B5^-2) ) * (( ( 1 - B5 * Ta ) * exp(-B5 * TD) ) - (1-B5*Tb)) +
    ( Cs * Vm * (B5^-2) ) * (( ( 1 - B5 * Ta ) * exp(-B5 * TD) ) - (1-B5*Tb)) +
    ( Cs * Vs * (B5^-3) ) * (( (-2 + (2*B5*Ta) - (B5^2*Ta^2)) * exp(-B5 * TD)) -
                             (-2 + (2*B5*Tb) - (B5^2 * Tb^2)))
