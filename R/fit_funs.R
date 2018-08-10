# Functions to fit a model and tune parameters by optimization


get_FEV1_prediction = function(df, Dos, K, A)

{
    dd = by(df, df$subject, function(x){
        browser
        experimentFEV1(x$O3, x$Ve, x$t_stop, Dos = Dos, K = K, A = A)})

    do.call(rbind, dd)
}

get_error = function(pred, obs, fun = MSE)
{
    fun(pred, obs)
}

MSE = function(x, y)
{
    (x-y)^2
}


fit_FEV1 = function(pars = c(Dos = 1100, K = 0.02, A = -0.02), df)
{
    pred = get_FEV1_prediction(df, Dos = pars[1], K = pars[2], A = pars[3])
    ans = sapply(seq(nrow(pred)), function(i){
        get_error(pred[i,], df$FEV1[df$subject == i])

    })
    mean(ans, na.rm = TRUE) #not quite right, but works for now 
}



