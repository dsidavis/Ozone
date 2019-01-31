## Tests that the model functions return exactly the same numbers
## as Ed's Excel sheet "Pred modV5.xlsx"
## Cols AU - BG

source("R/model_functions.R")

time = c(0,50, 60, 110, 120, 170, 215, 265, 275, 325, 335, 385, 399)
ve = rep(c(38.3, 15.0792453), 6)
o3 = c(0.0793933333, 0.0793933333, 0.07981, 0.07981,
       0.08029, 0.08029,
       0.07958, 0.07958,
       0.0803966666666, 0.0803966666666,
       0.07977333333333, 0.07977333333333) 

ans = experimentFEV1(o3, ve, time, Dos = 1400, A = -0.0246633, K = log(2)/35)

expected_ans = c(0, 0, 0, 0, 0, 0, -1.72087088354959e-94, -1.90017077526133e-49, 
-4.43485912124668, -4.19394040921599, -6.22098178707323, -5.47539345520911
)

stopifnot(all.equal(ans,expected_ans))

# O3 is too low to cause reaction
o3 = c(0.0005, 0.0005, 0.0004, 0.0004, 0.0005, 0.0005, 0.0005, 0.0005,
       0.0006, 0.0006, 0.0005, 0.0005)

ans2 = experimentFEV1(o3, ve, time, Dos = 1400, A = -0.0246633, K = log(2)/35)

stopifnot(all(ans2 == 0))
