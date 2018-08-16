# Functions to fit a model and tune parameters by optimization


get_FEV1_prediction = function(df, Dos, K, A)

{
    dd = by(df, list(df$person, df$protocolNum), function(x){
        # browser()
        experimentFEV1(x$O3, x$VE, x$endTime, Dos = Dos, K = K, A = A)})

    do.call(rbind, dd)
}

get_error = function(pred, obs, fun = MSE)
{
    fun(pred, obs)
}

SSE = function(x, y)
{
    (x-y)^2
}

loglike = function(x)
{
    # finish this later
    # dnorm()
}


fit_FEV1 = function(pars = c(Dos = 1100, K = 0.02, A = -0.02), df)
{
    pred = get_FEV1_prediction(df, Dos = pars[1], K = pars[2], A = pars[3])
    # browser()
    
    # ans = sapply(seq(nrow(pred)), function(i){
        # get_error(pred[i,], df$FEV1[df$person == i])

    # })
    
    sum(SSE(pred, df$dFEV1), na.rm = TRUE) #not quite right, but works for now 
}


fit_FEV1b = function(pars, d)

{
    sum(sapply(d, function(x, pars) fit_FEV1(pars, df = x), pars = pars), na.rm = TRUE)
}
