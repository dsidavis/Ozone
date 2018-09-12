ucd = list.files(path = "data", pattern = "UCD_.*\\.xls", full = TRUE)


fev = readSheet("data/UCD_WCA2002.xls")
o3 = readSheet("data/UCD_WCA2002.xls", "O3")
ans = combineFevO3(, fev, o3)


ans = combineFevO3("data/UCD_WCA2002.xls")

# Seems to work for
#   UCD_WCA2002.xls - manually checked 5 rows
#   UCD_WCA1997.xls - NAs for O3 correspond to the <i>-BS rows for all person values (1:12)
#   UCD_WCA2000b.xls - same as 1997
#
# Warnings for UCD_WCA2006ab.xls - getting dates/times when reading J50, P50 .. AI50 (21 cells).

# TO RESOLVE:
#    UCD_WCA2003ab.xls - NAs for O3 corresponding to -BS but also 1-1h, 1-2h, ... 3-r14
#                The O3 values are for 1-all, 2-all, 3-r, ...

#    UCD ESS2009 DATA.xls - different format.
