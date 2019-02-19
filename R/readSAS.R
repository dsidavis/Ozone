
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

mungeSAS_forStan = function(df)
    # Who is Stan? :P
{
    tmp = split(df, paste(df$ID, df$STUDY, df$Lab, df$EXPOSURE))
    # browser()

    t_vars = lapply(tmp, extract_t_vars)
    ind_vars = as.data.frame(do.call(rbind, lapply(tmp, extract_ind_vars)))
    dFEV1 = lapply(tmp, extract_dFEV1)

    t_vars = collapse_results(t_vars)
    n_timepts = lapply(t_vars, function(x) 
        apply(x, 2, function(x) sum(!is.na(x))))
    n_timepts = apply(do.call(rbind, n_timepts), 2, unique)
    t_vars = lapply(t_vars, trim_vars)

    dFEV1 = collapse_results(dFEV1)
    dFEV1$DELFEV1[is.na(dFEV1$DELFEV1)] = 0 #Stan cannot handle NA
    
    n_obs = length(tmp)
    n_ind = length(unique(ind_vars$ID))
    n_dFEV1 = sapply(tmp, nrow)
    
    list(max_timepts = max(n_timepts),
         max_n_dFEV1 = max(n_dFEV1),
         n_obs = n_obs,
         n_ind = n_ind,
         n_dFEV1 = n_dFEV1,
         n_timepts = n_timepts,
         ind = to_id(ind_vars$ID),
         age = ind_vars$AGE,
         BMI = ind_vars$BMI,
         BSA = ind_vars$BSA,
         Ve = t_vars$Ve,
         Cm = t_vars$Cm,
         Cs = t_vars$Cs,
         Time = t_vars$Time,
         dFEV1_measure_idx = dFEV1$TIME_ID,
         dFEV1 = dFEV1$DELFEV1)
    
}

to_id = function(x)
{
    as.integer(as.factor(x))
}



collapse_results = function(result_list)
    # takes a list of list, and collapses to list of df
{
    
    ans = lapply(seq_along(result_list[[1]]), function(i) {
        tmp = lapply(result_list, "[[", i)
        if(length(unique(sapply(tmp, length))) != 1)
            tmp = pad_list(tmp)
        do.call(cbind, tmp)    
    })
    names(ans) = names(result_list[[1]])
    ans
}

trim_vars = function(vars)
{
    i = apply(vars, 1, function(x) all(is.na(x)))
    vars[!i,]
}

pad_list = function(x, pad_value = 0)
{
    lens = sapply(x, length)
    max_len = max(lens)
    n_pad = max_len - lens
    lapply(seq_along(lens), function(i)
        c(x[[i]], rep(pad_value, n_pad[i])))
}

extract_t_vars = function(df,
                        varID = c("T", "Ve", "BSA",
                                  "O3_mean","O3_slope"),
                        varID_regex = paste0("^", varID, "_[0-9]{1,2}"))
{
    ans = lapply(varID_regex, function(regex)
        as.numeric(df[1, grep(regex, colnames(df))]))
    names(ans) = varID
    ans
}

extract_ind_vars = function(df,
                            varID = c("ID", "AGE", "BMI", "Male"))
{
    ans = as.numeric(df[1,varID])
    names(ans) = varID
    ans
}

extract_dFEV1 = function(df,
                         varID = c("DELFEV1", "TIME_ID"))
{
    df[, varID]
}

