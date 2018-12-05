x = c(1, 1, 1)
t = c(1, 3, 5)

#expect .5, 2.5, 4.5
cuml_integral(x, t)

#expect 1, 3, 5
cuml_integral(c(1, x), c(0, t))

# expect 0 0 0
cuml_integral(c(0,0,0), 1:3)

# triangle - expect last value == 3*3/2
cuml_integral(0:3, 0:3)
