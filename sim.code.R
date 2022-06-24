

install.packages("deSolve")
library(deSolve)
#return-get set up to save on github
#two species, interact only via density dependence, birth and death rate.

model <- function (time, y, parms) {
  with(as.list(c(y, parms)), {
    f <- P/(kP + P)
    dAlg <- r * f * Alg 
    dP   <- - r * 1/Y * f * Alg
    list(c(dAlg, dP))
  })
}
y      <- c(Alg = 0.1, P = 0.2)         # in mg/L
parms  <- c(r = 0.1, kP = 5e-3, Y = 41) # Y = C:P mass ratio
times  <- seq(0, 100, 1)
out <- ode(y, times, model, parms)
plot(out)
