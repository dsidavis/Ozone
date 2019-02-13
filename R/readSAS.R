
readSAS = function(file)
{
    d = as.data.frame(haven::read_sas(file))
    # browser()
    ans = by(d, list(d$ID, d$EXPOSURE),  mungeSAS)
    ans[!sapply(ans, is.null)]
}

mungeSAS = function(df)
{
    time = unlist(df[1,grep("^T_[0-9]{1,2}$", colnames(df))])
    ve = unlist(df[1,grep("^Ve_[0-9]{1,2}$", colnames(df))])
    o3 = unlist(df[1,grep("^O3_mean_[0-9]{1,2}$", colnames(df))])
    dfev1 = rep(NA, length(time))
    dfev1[df$TIME_ID -1] = if(all(is.na(df$DELFEV1))) 0 else df$DELFEV1 
    ans = data.frame(person = df$ID[1],
                     endTime = time,
                     VE = ve,
                     O3 = o3,
                     dFEV1 = -1 * dfev1/100, # need on decimal scale
                     sex = df$Male[1],
                     lab = df$Lab[1],
                     protocolNum = df$STUDY[1],
                     exposure = df$EXPOSURE[1],
                     stringsAsFactors = FALSE)
    rownames(ans) = NULL
    if(all(is.na(ans$endTime)))
        browser()
       

    ans[!is.na(ans$endTime),]
    
}
