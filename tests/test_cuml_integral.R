x = c(1, 1, 1)
t = c(1, 3, 5)

#expect .5, 2.5, 4.5
cuml_integral(x, t)

#expect 1
cuml_integral(c(1, x), c(0, t))
