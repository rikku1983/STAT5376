library(fda)
daytime=(1:365)-0.5
JJindex=c(182:365,1:181)
tempmat=daily$tempav[JJindex,]
tempbasis=create.fourier.basis(c(0,365), 65)
tempfd=smooth.basis(daytime,tempmat,tempbasis)$fd
tempfd$fdnames=list("Day (July 2 to June 30)",
                    "Weather Station",
                    "Mean temperature (deg. C)")
plot(tempfd, col=1,lty=1)#fig4.1


basis13 = create.bspline.basis(c(0,10), 13)
tvec = seq(0,1,len=13)
sinecoef = sin(2*pi*tvec)
sinefd = fd(sinecoef, basis13, list("t","","f(t)"))
op = par(cex=1.2)
plot(sinefd, lwd=2)
points(tvec*10, sinecoef, lwd=2)
par(op)

thawdata <- t(MontrealTemp[, 16:47])
daytime = ((16:47)+0.5)
par(cex=1.2)
plot(daytime, apply(thawdata,1,mean), "b", lwd=2,
     xlab="Day", ylab="Temperature (deg C)")
thawbasis=create.bspline.basis(c(16,48), nbasis = 7, norder=4)
thawbasismat=eval.basis(daytime,thawbasis)
thawcoef=solve(crossprod(thawbasismat),crossprod(thawbasismat,thawdata))
thawfd=fd(thawcoef,thawbasis,list('Day','Year','Temperature(C)'))
plot(thawfd,lty=1,lwd=2,col=1)#fig4.4
plotfit.fd(thawdata[,1],daytime,thawfd[1], lty=1,lwd=2)#4.5

#Linear Differential Operator: Lfd class
betalist=vector('list',3)
thawconst.basis=create.bspline.basis(c(16,48),nbasis=1,norder=1)
betalist[[1]] = fd(0, thawconst.basis)
betalist[[2]] = fd(2, thawconst.basis)
betalist[[3]] = fd(0, thawconst.basis)
harmaccelLfd=Lfd(3,betalist)

harmaccelLfd = vec2Lfd(c(0,2,0), c(0, 365))
Ltempmat = eval.fd(daytime, tempfd, harmaccelLfd)

#Plot derivative
D1tempfd = deriv.fd(tempfd, 1)
plot(D1tempfd)
#Plot Linear differential operator
w=1
harmaccelLfd = vec2Lfd(c(0,w^2,0), c(0, 365))
Ltempfd = deriv.fd(tempfd, harmaccelLfd)
plot(Ltempfd)

#Bivariate Functional data objects
heightbasis12=create.bspline.basis(c(1,18),12,6)
rndbspl=create.bspline.basis(c(0,1), nbasis=23,norder = 4)
plot(rndbspl)
coef=(1:23)+rnorm(23)
rndbsplfd=fd(coef, rndbspl)
plot(rndbsplfd)
plotfit.fd(coef[c(1,3:21,23)], seq(0,1,0.05), rndbsplfd)
D1=deriv.fd(rndbsplfd,1)
D2=deriv.fd(rndbsplfd,2)
plot(c(rndbsplfd,D1,D2))
plot(D1)
plot(D2)

########################################################
##########Chapter 5 Computing Curves from Noisy Data############
###########################################################
heightmat=growth$hgtf
age=growth$age
heightbasis12=create.bspline.basis(c(1,18),12,6)
#age <- c( seq(1, 2, 0.25), seq(3, 8, 1), seq(8.5, 18, 0.5))
basismat=eval.basis(age,heightbasis12)
heightcoef=lsfit(basismat, heightmat, intercept = F)$coef
heightfd1=fd(heightcoef, heightbasis12, list('age','Girl','Height'))
# Off the shelf function in r
heightList=smooth.basis(age, heightmat, heightbasis12)
heightfd = heightList$fd
height.df = heightList$df
height.gcv = heightList$gcv
# Above methods are both regression splines
# the instability of regression spline derivative estimates at the
# boundaries is especially acute

#Roughness Penalties
# uses a large number of basis functions, possibly
# extending to one basis function per observation and even beyond
# D2 and harmonic acceleration

Rmat=eval.penalty(tempbasis,harmaccelLfd)
# As l !0;df(l )!min(n;K), where n = the number of observations and K =
#   the number of basis functions. Similarly, as l ! ¥;df(l ) ! m, where m is the
# order of the highest derivative used to define the roughness penalty.

# Functional Parameter Objects
norder=6
nbasis=length(age)+norder-2
heightbasis=create.bspline.basis(c(1,18), nbasis, norder, age)
heightfdPar=fdPar(heightbasis, Lfdobj=4, lambda=0.01)
heightfd=smooth.basis(age,heightmat,heightfdPar)$fd
plot(heightfd)

#5.2.5 choosing lambda
# by GCV
loglam=seq(-6,0,0.25)
gcvsave=matrix(0, nrow=length(loglam), ncol=1)
dfsave=gcvsave
for(i in 1:length(loglam)){
  lambdai=10^loglam[i]
  hgtfdPari=fdPar(heightbasis, 4, lambdai)
  smthfd=smooth.basis(age, heightmat, hgtfdPari)
  gcvsave[i]=sum(smthfd$gcv)
  dfsave[i]=smthfd$df
}
plot(loglam, gcvsave, 'o', lwd=2)

#Case study
dayOfYearShifted=c(182:365,1:181)
logprecav=CanadianWeather$dailyAv[dayOfYearShifted,,'log10precip']
#set up a saturated Fourier basis for the data
dayrange=c(0,365)
daybasis = create.fourier.basis(dayrange, 365)
Lcoef=c(0,(2*pi/diff(dayrange))^2,0)
harmaccelLfd=vec2Lfd(Lcoef, dayrange)
loglam=seq(4,9,0.25)
nlam=length(loglam)
dfsave=rep(NA, nlam)
gcvsave=rep(NA,nlam)
for(ilam in 1:nlam){
  cat(paste('log10 lambda =',loglam[ilam],'\n'))
  lambda=10^loglam[ilam]
  fdParobj=fdPar(daybasis, harmaccelLfd, lambda)
  smoothlist=smooth.basis(day.5,logprecav, fdParobj)
  dfsave[ilam]=smoothlist$df
  gcvsave[ilam]=sum(smoothlist$gcv)
}
plot(loglam, gcvsave, 'o', lwd=2)
#Pick lambda at 10 to 6 power which minimize the gcv
lambda=1e6
fdParobj=fdPar(daybasis, harmaccelLfd, lambda)
logprec.fit=smooth.basis(day.5, logprecav, fdParobj)
logprec.fd=logprec.fit$fd
fdnames=list('Day(7/1 to 6/30)','Weather Station'=CanadianWeather$place,
             'Log 10 precipitation (mm)')
logprec.fd$fdnames=fdnames
plot(logprec.fd)
plotfit.fd(logprecav[3], day.5, logprec.fd[3])

#Other constraints, by using transformation

# Assessing Fit
logprecmat=eval.fd(day.5, logprec.fd)
logprecres=logprecav-logprecmat
# across station
logprecvar1=apply(logprecres^2,1,sum)/35
logprecvar2=apply(logprecres^2,2,sum)/(365-12)

plot(logprecvar2^(0.5)) #fig5.7
logstddev.fd=smooth.basis(day.5, log(logprecvar1)/2, fdParobj)$fd
logprecvar1fit=exp(eval.fd(day.5,logstddev.fd))
plot(day.5, logprecvar1fit,col='red','l',ylim=c(0,0.25))
points(day.5, logprecvar1^0.5)#fig5.8


#CH6 Description of functions
plot(logprec.fd)
meanlogprec=mean(logprec.fd)
plot(meanlogprec, lwd=4,add=T)
stddevlogprec=std.fd(logprec.fd)
str(stddevlogprec)
plot(stddevlogprec)
logprecvar.bifd = var.fd(logprec.fd)
weektime = seq(0,365,length=53)
logprecvar_mat = eval.bifd(weektime, weektime,
                           logprecvar.bifd)
persp(weektime, weektime, logprecvar_mat, theta = -45,
      phi=25,r=3,expand=0.5) #fig6.1
contour(weektime,weektime, logprecvar_mat) #fig6.2

#Probes like contrast in Experimental design
# probeval = inprod(xifd, xfd)

# Phase-plane plots of periodic effects
agefine=seq(1,18,length.out=101)
velffine=eval.fd(agefine, heightfd[1:10],1)
accffine=eval.fd(agefine,heightfd[1:10],2)
plot(velffine[,1:5],accffine[,1:5],'o')
points(velffine[64,1:5],accffine[64,1:5],col='red') #fig1.15

#Confidence Interval(CI)
dayvec=seq(0,365,len=101)
xivec=exp(20*cos(2*pi*(dayvec-197)/365))
xibasis=create.bspline.basis(c(0,365),13)
xifd=smooth.basis(dayvec,xivec,xibasis)$fd
tempbasis=create.fourier.basis(c(0,365), 65)
precbasis=create.fourier.basis(c(0,365), 365)
tempLmat=inprod(tempbasis,xifd)
precLmat=inprod(precbasis,xifd)

#Plot CL for station 'Prince Rupert'
lambda=1e6
fdParobj=fdPar(daybasis, harmaccelLfd, lambda)
logprecList=smooth.basis(day.5,logprecav,fdParobj)
logprec.fd=logprecList$fd
fdnames=list("Day (July 1 to June 30)",
             "Weather Station" = CanadianWeather$place,
             "Log 10 Precipitation (mm)")
logprec.fd$fdnames=fdnames
#Estimate Covar matrix of residuals
# Next we estimate Se, which we assume is diagonal. Consequently, we need only
# estimate the variance of the residuals across weather stations for each day. We do
# this by smoothing the log of the mean square residuals and then exponentiating the
# result:
logprecmat=eval.fd(day.5,logprec.fd)
logprecres=logprecav-logprecmat
logprecvar=apply(logprecres^2,1,sum)/(35-1)
lambda=1e8
resfdParobj=fdPar(daybasis,harmaccelLfd, lambda)
logvar.fit=smooth.basis(day.5, log(logprecvar), resfdParobj)
logvar.fd=logvar.fit$fd
varvec=exp(eval.fd(day.5,logvar.fd))
SigmaE=diag(as.vector(varvec))

# Next we get y2cMap from the output of smooth.basis, and compute c2rMap
# by evaluating the smoothing basis at the sampling points. We then compute the
# variance-covariance matrix for curve values, and finish by plotting the log precipitation
# curve for Prince Rupert along with this curve plus and minus two standard
# errors. The result is Figure 6.6.
y2cMap=logprecList$y2cMap 
c2rMap=eval.basis(day.5, daybasis)
Sigmayhat=c2rMap %*% y2cMap %*% SigmaE %*% t(y2cMap) %*% t(c2rMap)
logprec.stderr=sqrt(diag(Sigmayhat))
logprec29=eval.fd(day.5,logprec.fd[29])
plot(logprec.fd[29],lwd=2,ylim=c(0.2,1.3))
lines(day.5,logprec29+2*logprec.stderr,lty=2,lwd=2)
lines(day.5,logprec29-2*logprec.stderr,lty=2,lwd=2)
points(day.5, logprecav[,29]) #fig6.6
