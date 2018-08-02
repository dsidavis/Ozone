library(readxl)

get_times = function(xls)
{
    i = grep("^TIME$", excel_sheets(xls))
    if(length(i) == 0)
        stop("Time sheet not found")

    if(length(i) > 1)
        stop("More than one 'Time' sheet found")
    
    as.data.frame(read_excel(xls, sheet = i), stringsAsFactors = FALSE)
}
