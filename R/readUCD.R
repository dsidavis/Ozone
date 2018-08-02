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
