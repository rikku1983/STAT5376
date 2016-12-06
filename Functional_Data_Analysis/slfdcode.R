library(fda)
#background
age = growth$age
#  define the range of the ages and set up a fine mesh of ages
age.rng = range(age)
agefine = seq(age.rng[1], age.rng[2], length=501)
# Monotone smooth (see section 5.4.2)
# B-spline of order 6 = quintic polynomials
# so the acceleration will be cubic
gr.basis = create.bspline.basis(norder=6, breaks=growth$age)
# Consider only the first 10 girls
children = 1:10
ncasef   = length(children)
# matrix of starting values for coefficients
cvecf           = matrix(0, gr.basis$nbasis, ncasef)
dimnames(cvecf) = list(gr.basis$names,
                       dimnames(growth$hgtf)[[2]][children])

# Create an initial functional data object with these coefficients

gr.fd0  = fd(cvecf, gr.basis)

# Create an initial functional parameter object
# Lfdobj = 3 to penalize the rate of change of acceleration

gr.Lfd    = 3

# Figure 1.15 was created with lambda = 10^(-1.5);
# we use that also for Figure 1.1

gr.lambda = 10^(-1.5)

#  Define the functional parameter object

gr.fdPar  = fdPar(gr.fd0, Lfdobj=gr.Lfd, lambda=gr.lambda)

#  Monotonically smooth the female data

hgtfmonfd   = with(growth, smooth.monotone(age, hgtf[,children], gr.fdPar))

#  Plot Figure 1.1

#   The with function evaluates an R expression in an environment constructed
#   from data, possibly modifying the original data.

with(growth, matplot(age, hgtf[, children], pch='o', col=1,
                     xlab='Age (years)', ylab='Height (cm)',
                     ylim=c(60, 183)) )

hgtf.vec = predict(hgtfmonfd$yhatfd, agefine)
matlines(agefine, hgtf.vec, lty=1, col=1)