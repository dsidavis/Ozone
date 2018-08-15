readUCD = 
function(xls = "UCD_WCA2002.xls", data = readxl::read_excel(xls, "FEV1"))
{
    d = as.data.frame(data)
    fev = unlist(d[-(1:5), -(1:5)])

    vars = gsub(" ", "", d[5, 1:5])
    d = d[-(1:5), 1:5]
    names(d) = vars
    
    numPeople = ncol(data) - 5L

    data.frame(fev1 = fev,
               startTime = rep(d$StartTime, numPeople),
               endTime = rep(d$EndTime, numPeople),
               period = rep(d$SamplePeriod, numPeople),
               protocolNum = rep(d$"Protocol#", numPeople),
               Descriptor = rep(d$Descriptor, numPeople),
               person = rep(seq_len(numPeople), each = nrow(d)),
               study = rep(gsub("\\.xls$", "", xls), length(fev)))
}


readSheet =
function(xls = "UCD_WCA2002.xls", sheet = "FEV1", data = readxl::read_excel(xls, sheet))
{
    d = as.data.frame(data)
    fev = unlist(d[-(1:5), -(1:5)])

    vars = gsub(" ", "", d[5, 1:5])
    d = d[-(1:5), 1:5]
    names(d) = vars
    
    numPeople = ncol(data) - 5L

    ans = data.frame(fev1 = fev,
               startTime = rep(d$StartTime, numPeople),
               endTime = rep(d$EndTime, numPeople),
               period = rep(d$SamplePeriod, numPeople),
               protocolNum = rep(d$"Protocol#", numPeople),
               # Order the descriptors as 1-BS, 1-1h, 1-2h, ... and then 2-BS, ...
               # Make an ordered factor.
               Descriptor = rep(gsub(" -", "-", as.character(d$Descriptor)), numPeople),
               person = rep(seq_len(numPeople), each = nrow(d)),
               study = rep(gsub("\\.xls$", "", xls), length(fev)), stringsAsFactors = FALSE)

    names(ans)[1] = sheet
    
    invisible(ans)
}

readExp = function(xls = "UCD_WCA2002.xls", sheets = c("FEV1", "EXVE", "O3", "RVE"),
                   data = lapply(sheets, function(x, xls) readSheet(xls, sheet = x), xls = xls))
{
    d = lapply(data, as.data.frame)

    # assumes all sub. accounted for
    n_person = sapply(d, function(x) length(unique(x$person))) 
    if(any(diff(n_person) != 0))
        stop("Different number of subjects between sheets!\n",
             paste(sheets, n_person, sep = ": ", collapse = "\n"))
    browser()

    n_protocol = sapply(d, function(x) length(unique(x$protocolNum)))

    if(any(diff(n_protocol) != 0))
        stop("Different number of protocol between sheets!\n",
             paste(sheets, n_protocol, sep = ": ", collapse = "\n"))

    ans = lapply(seq(n_person[1]), function(i) {
        lapply(seq(n_protocol[1]), function(j) {
            dd = lapply(d, function(x) {
                tmp = x[x$person == i & x$protocolNum == j,]
                expandTS(tmp$startTime, tmp$endTime, tmp[,1], colnames(tmp)[1])
            })
            a = Reduce(function(x, y) merge(x, y, all = TRUE), dd)
            a$person = i
            a$protocolNum = j
            a
        })
    })
    
    ans = unlist(ans, recursive=FALSE)
    ans = lapply(ans, summarizeTS)
    
    d
    
}

expandTS = function(starts, ends, value, val_name)
{
    
    starts = as.integer(starts)
    ends = as.integer(ends)

    ans = lapply(seq_along(starts), function(i) {
        
        time = starts[i]:max(c(ends[i]-1, 0)) #avoid overlapping times
        data.frame(time = time,
                   value = value[i])
    })
    ans = do.call(rbind, ans)
    colnames(ans)[2] = val_name
    ans
}

summarizeTS = function(ts)
{
    ts = combineVE(ts)
    i = !duplicated(ts[,colnames(ts) != "time"])
    ts_reduced = ts[i,]

    ans = lapply(seq(nrow(ts_reduced)), function(j){
        tmp = ts_reduced[j,]
        tmp$startTime = tmp$time
        if(j < nrow(ts_reduced)){
            tmp$endTime = ts_reduced$time[j+1]
        } else
            tmp$endTime = max(ts_reduced$time)
        tmp
    })
    do.call(rbind, ans)
}

combineVE = function(df, dropExtra = TRUE)
{
    exve = !is.na(df$EXVE)
    rve = !is.na(df$RVE)
    if(any(exve & rve))
        warning("Values present for both EXVE and RVE. Using RVE")
    if(any(!exve & !rve))
        warning("Value missing for both EXVE and RVE. Using NA")
    ve = df$EXVE
    ve[!exve] = df$RVE[rve]
    df$VE = ve
    if(dropExtra)
        df = df[, !colnames(df) %in% c("EXVE", "RVE")]
    
    df
}


mkId =
function(d, vars)
{
  apply(d[, vars], 1, paste, collapse = ",")
}




combineFevO3 = 
function(xls = "UCD_WCA2002.xls", fev = readSheet(xls), o3 = readSheet(xls, "O3"),
         addOzoneTimes = TRUE, addBProtocol = FALSE)
{        

   i = match(mkId(fev, c("Descriptor", "person")), mkId(o3, c("Descriptor", "person")), 0)
   fev$O3 = NA
   fev[ i != 0, "O3" ] = o3$O3[i]

#fev[is.na(fev$O3),]

    # Do we need the start and endTime from 03
    if(addOzoneTimes) {
         # For now add personO3 so we can verify the matched.
        fev$DescriptorO3 = fev$personO3 = fev$endTimeO3 = fev$startTimeO3 = NA
        fev[ i != 0, c("startTimeO3", "endTimeO3", "personO3", "DescriptorO3") ] = o3[i, c("startTime", "endTime", "person", "Descriptor") ]
    }

# So the <i>-BS are not in the O3 data and we get O3= NA for these
# and we don't have the <i>-<j>b rows from the O3 in fev since there is
    # If we wanted to add the <i>-<j>b rows tpo fev, we can do the following:

    if(addBProtocol) {
        j = match(o3$Desc, fev$Desc)
        # table(droplevels(o3[is.na(j), "Descriptor"]))

        tmp = o3[is.na(j),]
        tmp2 = cbind(FEV1 = rep(NA, nrow(tmp)), tmp)
        tmp = rbind(fev, tmp2[names(fev)])
    } else
        fev

}    

