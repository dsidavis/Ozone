x = c(1, 1, 1)
t = c(1, 3, 5)

#expect .5, 2.5, 4.5
stopifnot(identical(cuml_integral(x, t), c(.5,2.5, 4.5)))

stopifnot(identical(cuml_integral(x*2, t), c(.5,2.5, 4.5) * 2))

#expect 1, 3, 5
stopifnot(identical(cuml_integral(c(1, x), c(0, t)), c(0,1, 3, 5)))

# expect 0 0 0
stopifnot(identical(cuml_integral(c(0,0,0), 1:3), rep(0,3)))

# triangle - expect last value == 3*3/2
stopifnot(cuml_integral(0:3, 0:3)[4] == 3*3/2) 
