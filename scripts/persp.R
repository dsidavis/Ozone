# The idea here is to look at the sum of squares that we computed from the grid
# of parameter values in gridSearch.R

# We read the R object from the RDS file. This is a data frame
# with 100^3 rows corresponding  to a grid of 100 unique values for each
# of the 3 variables DOS, K and A. The FEV1 value is the Sum of Squared errors.

# rds = "~/VeryFineGrained_200.rds"
rds = "grid.rds"
                                        #rds = "FinerGrain.rds"
p = readRDS(rds)

# keep only the rows that are within 5% of the minumum
w = p$SSE < min(p$SSE)*1.1
# .8% of the values are within 5%
pm = p[w,]

plot(density(pm$SSE), xlim = range(pm$SSE))

# Considering the range of values for SSE (FEV1) (2.11 to 13 million)
# the range we see within 5% of the minimum is very small. So 8000 values within .1 of the minimum.


plot(density(pm$DOS))
plot(density(pm$DOS)); rug(pm$DOS)
plot(SSE ~ DOS, pm)


# Turn the pm data frame into a 3-dimensional array of SSE values. This is for use with persp().
# We create the 3-D array with dimnames corresponding to the unique values of DOS, K and A in pm.
# Then we map the rows in pm to the cells in the array.
# Compute the dimnames from the unique values of each of the three parameters, to 8 digits.

u = lapply(pm[1:3], function(x) formatC(sort(unique(x)), 8))
a = array(as.numeric(NA), sapply(u, length), dimnames = u)

# Create a 3 column matrix giving the names of the indices in a so that we can use
# this to set the FEV1 value for each row in pm into the corresponding cell in a.
i = do.call(cbind, lapply(pm[1:3], formatC, 8))

a[i] = pm$SSE

par(mfrow = c(3, 5))
invisible(sapply(1:dim(a)[3], function(i) persp(a[,,i], xlab = "A", main = dimnames(a)[[3]][i])))




plot(SSE ~ DOS, pm, col = heat.colors(10)[ cut(K, 10) ])


ggplot(pm) + geom_point(aes(x = DOS, y = SSE, color = K)) + scale_color_viridis_c()

ggplot(pm) + geom_point(aes(x = DOS, y = SSE, color = A))

ggplot(pm) + geom_point(aes(x = DOS, y = SSE, color = K)) + scale_color_viridis_c()
ggplot(pm) + geom_point(aes(x = K, y = SSE, color = DOS)) + scale_color_viridis_c()
ggplot(pm) + geom_point(aes(x = K, y = A, color = SSE)) + scale_color_viridis_c()
ggplot(pm) + geom_point(aes(x = K, y = DOS, color = SSE)) + scale_color_viridis_c()


library(plot3D)
#with(pm, scatter3D(DOS, K, A, bty = "g", pch = 12, col = gg.col(100)))
with(pm, scatter3D(DOS, K, A, bty = "g", pch = 12, colvar = SSE, xlab = "DOS", ylab = "K", zlab = "A"))

with(pm, scatter3D(DOS, K, A, bty = "g", pch = 12, colvar = SSE, xlab = "DOS", ylab = "K", zlab = "A", phi = 10, theta = 10))
