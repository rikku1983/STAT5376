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
thawbasismat=eval.dimbasis(daytime,thawbasis)
thawcoef=solve(crossprod(thawbasismat),crossprod(thawbasismat,thawdata))
thawfd=fd(thawcoef,thawbasis,list('Day','Year','Temperature(C)'))
plot(thawfd,lty=1,lwd=2)#fig4.4
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
w=2*pi/32
harmaccelLfd = vec2Lfd(c(0,w^2,0), c(0, 365))
Ltempfd = deriv.fd(thawfd, harmaccelLfd)
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

#Understanding Ldf
x=seq(0,2*pi,len=200)
y=sin(x)
plot(x,y)
sinbasis=create.fourier.basis(c(0,2*pi), 21)
slsin=smooth.basis(x,y,sinbasis)$fd
plot(slsin)
sllfd1=vec2Lfd(c(0,10),c(0,2*pi))
sllfd2=vec2Lfd(c(0,5),c(0,2*pi))
sllfd3=vec2Lfd(c(0,1,1),c(0,2*pi))
slvel1=eval.fd(x,slsin,sllfd1)
slvel2=eval.fd(x,slsin,sllfd2)
plot(x,slvel1,'l')
lines(x,slvel2,'l',col='red')
slacc=eval.fd(x,slsin,harmaccelLfd)
plot(slvel,slacc)

#CH7 FPCA
plot(logprecList$fd[2]-mean(logprecList$fd))
mlogprec<-mean(logprecList$fd)
class(mlogprec)
plot(mlogprec)

logprecfd=logprecList$fd
logprec.pcalist=pca.fd(logprecfd,2)
print(logprec.pcalist$values)
plot.pca.fd(logprec.pcalist)
logprec.rotapcalist=varmx.pca.fd(logprec.pcalist)
plot.pca.fd(logprec.rotapcalist)

#PCA on residuals
logprecres.fd=smooth.basis(day.5, logprecres, fdParobj)$fd
plot(logprecres.fd, lwd=2, col=1, lty=1, cex=1.2,
     xlim=c(0,365),ylim=c(-0.07,0.07), xlab='Day',
     ylab='Residual')#fig 7.4
logprecres.pcalist=pca.fd(logprecres.fd, 2)
plot.pca.fd(logprecres.pcalist)
#More PCR
heightfd1pca=pca.fd(heightfd1,2)
heightfd1pca$harmonics$coefs
heightfd1$coefs
heightfd1pca$values
plot.pca.fd(heightfd1pca)
heightfd1rotpca=varmx.pca.fd(heightfd1pca)
plot.pca.fd(heightfd1pca)
heightm=eval.fd(age,heightfd1)
heightres=heightmat-heightm
heightresfd=smooth.basis(age, heightres,heightbasis12)$fd
plot(heightresfd)
heightrespca=pca.fd(heightresfd,2)
plot.pca.fd(heightrespca)
heightresrotpca=varmx.pca.fd(heightrespca)
plot.pca.fd(heightrespca)

#PCA on handwritten fda
dim(handwrit)
fdarange=c(0,2300)
fdabasis=create.bspline.basis(fdarange,105,6)
fdatime=seq(0,2300,len=1401)
fdafd=smooth.basis(fdatime,handwrit, fdabasis)$fd
fdafd$fdnames[[1]]='Milliseconds'
fdafd$fdnames[[2]]='Replications'
fdafd$fdnames[[3]]=list('x','y')
plot(fdafd)
#PCA
nharm=3
fdapcaList=pca.fd(fdafd,nharm)
plot.pca.fd(fdapcaList)
fdarotpcaList=varmx.pca.fd(fdapcaList)
plot.pca.fd(fdarotpcaList)

fdaeig=fdapcaList$values
neig=12
x=matrix(1,neig-nharm,2)
x[,2]=(nharm+1):neig
y=log10(fdaeig[(nharm+1):neig])
c=lsfit(x,y,int=F)$coef
par(mfrow=c(1,1),cex=1.2)
plot(1:neig, log10(fdaeig[1:neig]),'b',
     xlab='Eigenvalue Number',
     ylab='Log10 Eigenvalue')
lines(1:neig,c[1]+c[2]*(1:neig),lty=2)

#Fig7.7
fdameanfd  = mean(fdafd)
fdameanmat = eval.fd(fdatime, fdameanfd)
#  evaluate the harmonics
harmfd  = fdarotpcaList$harm
harmmat = eval.fd(fdatime, harmfd)

fdapointtime = seq(0,2300,len=201)
fdameanpoint = eval.fd(fdapointtime, fdameanfd)
harmpointmat = eval.fd(fdapointtime, harmfd)

fac = 0.1
harmplusmat = array(0,c(201,3,2))
harmminsmat = array(0,c(201,3,2))
for (j in 1:3) {
  harmplusmat[,j,] = fdameanpoint[,1,] + fac*harmpointmat[,j,]
  harmminsmat[,j,] = fdameanpoint[,1,] - fac*harmpointmat[,j,]
}
j=3
plot(fdameanmat[,1,1]-0.035,  fdameanmat[,1,2], "l", lwd=2,
     xlim=c(-0.075,0.075), ylim=c(-0.04, 0.04),
     xlab="", ylab="")
lines(harmplusmat[,j,1]-0.035, harmplusmat[,j,2], lty=2)
lines(harmminsmat[,j,1]-0.035, harmminsmat[,j,2], lty=2)
j=2
lines(fdameanmat[,1,1]+0.035,  fdameanmat[,1,2],  lty=1, lwd=2)
lines(harmplusmat[,j,1]+0.035, harmplusmat[,j,2], lty=2)
lines(harmminsmat[,j,1]+0.035, harmminsmat[,j,2], lty=2)

#CCA
tempav = CanadianWeather$dailyAv[
  dayOfYearShifted, , 'Temperature.C']
lambda   = 1e2
fdParobj = fdPar(daybasis, harmaccelLfd, lambda)
temp.fd  = smooth.basis(day.5, tempav, fdParobj)$fd
temp.fd$fdnames = list("Day (July 2 to June 30)",
                       "Weather Station",
                       "Mean temperature (deg. C)")
ccafdPar = fdPar(daybasis, 2, 5e6)
ccalist  = cca.fd(temp.fd, logprec.fd, 3, ccafdPar, ccafdPar)

ccawt.temp    = ccalist$ccawtfd1
ccawt.logprec = ccalist$ccawtfd2
corrs         = ccalist$ccacorr
print(corrs[1:3])
#  [1] 0.9139817 0.6194850 0.3495515
ccawtmat.temp    = eval.fd(day.5, ccawt.temp)
ccawtmat.logprec = eval.fd(day.5, ccawt.logprec)
#  Figure 7.8
plot(day.5, ccawtmat.temp[,1], type='l', lwd=2, cex=2,
     xlab="Day (July 1 to June 30)",
     ylab="Canonical Weight Functions")
lines(day.5, ccawtmat.logprec[,1], lty=2, lwd=2)
lines(dayrange, c(0, 0), lty=3)
legend("bottomleft", c("Temp.", "Log Prec."), lty=c(1,2))
#  Figure 7.9
ccascr.temp    = ccalist$ccavar1
ccascr.logprec = ccalist$ccavar2

placeindex = c(35,30,31,19,33,25,24,17,16,8,14,12,15,10,27,6,1,29)

plot(ccascr.temp[,1], ccascr.logprec[,1], type="p", pch="*", cex=2,
     xlim=c(-40,80),
     xlab="Temperature Canonical Weight",
     ylab="Log Precipitation Canonical Weight")
text(ccascr.temp[placeindex,1]+10, ccascr.logprec[placeindex,1],
     CanadianWeather$place[placeindex])


#CH8 Registration
#Growth data
##  Compute the monotone smoothing of the Berkeley female growth data.
##

#  set up ages of measurement and an age mesh

age     = growth$age
nage    = length(age)
ageRng  = range(age)
nfine   = 101
agefine = seq(ageRng[1], ageRng[2], length=nfine)

#  the data

hgtf   = growth$hgtf
ncasef = dim(hgtf)[2]

#  an order 6 bspline basis with knots at ages of measurement

norder = 6
nbasis = nage + norder - 2
wbasis = create.bspline.basis(ageRng, nbasis, norder, age)

#  define the roughness penalty for function W

Lfdobj    = 3          #  penalize curvature of acceleration
lambda    = 10^(-0.5)  #  smoothing parameter
cvecf     = matrix(0, nbasis, ncasef)
Wfd0      = fd(cvecf, wbasis)
growfdPar = fdPar(Wfd0, Lfdobj, lambda)

#  monotone smoothing

growthMon = smooth.monotone(age, hgtf, growfdPar)

# (wait for an iterative fit to each of 54 girls)

Wfd        = growthMon$Wfd
betaf      = growthMon$beta
hgtfhatfd  = growthMon$yhatfd

plot(hgtfhatfd)

accelfdUN=deriv.fd(hgtfhatfd,2)
plot(accelfdUN)
accelmeanfdUN=mean(accelfdUN)
PGSctr=rep(0,10)
agefine=seq(1,18,len=101)
par(mfrow=c(1,1),ask=T)
for(icase in 1:11){
  accveci=predict(accelfdUN[icase],agefine)
  plot(agefine,accveci,'l',ylim=c(-6,4),xlab="Year", ylab="Height Accel.",
       main=paste("Case",icase))
  lines(c(1,18),c(0,0),lty=2)
  PGSctr[icase]=locator(1)$x
}
PGSctrmean=mean(PGSctr)
wbasisLM=create.bspline.basis(c(1,18),4,3,c(1,PGSctrmean,18))
WfdLM=fd(matrix(0,4,1),wbasisLM)
WfdParLM=fdPar(WfdLM, 1, 1e-12)

regListLM=landmarkreg(accelfdUN[1:11], PGSctr, PGSctrmean, WfdParLM, TRUE)
accelfdLM     = regListLM$regfd
warpLM=regListLM$warpfd
accelmeanfdLM = mean(accelfdLM)
plot(accelfdLM)

#Countinuous registration
wbasisCR=create.bspline.basis(c(1,18),15,5)
Wfd0CR=fd(matrix(0,15,11), wbasisCR)
WfdParCR=fdPar(Wfd0CR,1,1)
regList=register.fd(mean(accelfdLM),accelfdLM, WfdParCR)
accelfdCR=regList$regfd
warpfdCR  = regList$warpfd
plot(accelfdCR)
#amp and phase decomposition
AmpPhasList=AmpPhaseDecomp(accelfdUN[1:11],accelfdLM, warpLM, c(3,18))
  


#CH9 Functional regression with scalar response
#predict annual precipitation using temperature profile
annualprec=log10(apply(daily$precav,2,sum))
#define x 35 temperature curves
tempbasis=create.fourier.basis(c(0,365),65)
tempSmooth=smooth.basis(day.5, daily$tempav,tempbasis)
tempfd=tempSmooth$fd
#define xfdlist
templist=vector('list',2)
templist[[1]]=rep(1,35)
templist[[2]]=tempfd

#9.4.1 Low-dimension regression coefficient function beta
conbasis=create.constant.basis(c(0,365))
betabasis=create.fourier.basis(c(0,365),5)
betalist=vector("list",2)
betalist[[1]]=conbasis
betalist[[2]]=betabasis
fRegressList=fRegress(annualprec,templist,betalist)
#Plot estimate of regression function
betaestlist=fRegressList$betaestlist
tempbetafd=betaestlist[[2]]$fd
plot(tempbetafd, xlab='day', ylab='beta for temp')
#assess the fit
annualprechat1=fRegressList$yhatfdobj
annualprecres1=annualprec-annualprechat1
SSE1.1=sum(annualprecres1^2) #SSE of our model
SSE0=sum((annualprec-mean(annualprec))^2)
RSQ1=(SSE0-SSE1.1)/SSE0
Fratio1=((SSE0-SSE1.1)/5)/(SSE1.1/29)

#Coef beta estimate using roughness penalty
#Lcoef=c(0,1,0)
#harmaccelLfd=vec2Lfd(Lcoef,c(0,365))
#test harmonic acceleration operator
#testbasis=create.fourier.basis(c(0,365),3)
#testfd=fd(c(1,1,1), testbasis)
# eval.basis(0,testbasis)
# eval.fd(0,testfd)
# eval.basis(1/2/2,testbasis)
# plot(testbasis)
# c1=sqrt(omega/(2*pi)) to make sure orthonormal functions
# testderiv=deriv.fd(testfd, harmaccelLfd0)
# plot(testderiv)
# 
# 
# hmat <- vec2Lfd(matrix(Lcoef, 1), c(0, 365))
# 
# 
# all.equal(harmaccelLfd, hmat)
# 
# omega = 2*pi/365
# thawconst.basis = create.constant.basis(c(0,365))
# 
# betalist = vector("list", 3)
# betalist[[1]] = fd(0, thawconst.basis)
# betalist[[2]] = fd(omega^2, thawconst.basis)
# betalist[[3]] = fd(0, thawconst.basis)
# #betalist[[4]] = fd(1, thawconst.basis)
# harmaccelLfd0 = Lfd(3, betalist)
# class(betalist[[3]])
Lcoef <- c(0,(2*pi/365)^2,0)
harmaccelLfd <- vec2Lfd(Lcoef, c(0,365))
betabasis=create.fourier.basis(c(0,365),35)
lambda=10^12.5
betafdPar=fdPar(betabasis, harmaccelLfd, lambda)
betalist[[2]]=betafdPar

annPrecTemp=fRegress(annualprec, templist, betalist)
betaestlist2=annPrecTemp$betaestlist
annualprechat2=annPrecTemp$yhatfdobj
print(annPrecTemp$df)
SSE1.2=sum((annualprec-annualprechat2)^2)
RSQ2=(SSE0-SSE1.2)/SSE0
Fratio2=((SSE0-SSE1.2)/3.7)/(SSE1.2/30.3)
#fig9.2
plot(annualprechat2, annualprec, lwd=2)
abline(lm(annualprec~annualprechat2), lty='dashed', lwd=2)

#Constant beta modal
# betalist[[2]]=fdPar(conbasis)
# fRegressList0=fRegress(annualprec,templist,betalist)
# betaestlist=fRegressList$betaestlist
# annualprechat = fRegressList$yhatfdobj
# SSE1 = sum((annualprec-annualprechat)^2)
# RSQ = (SSE0 - SSE1)/SSE0
# Fratio = ((SSE0-SSE1)/1)/(SSE1/33)

# 9.4.3 Choosing Smoothing Parameters

loglam = seq(5,15,0.5)
nlam   = length(loglam)
SSE.CV = rep(NA,nlam)
for (ilam in 1:nlam) {
  #print(paste("log lambda =", loglam[ilam]))
  lambda     = 10^(loglam[ilam])
  betalisti  = betalist
  betafdPar2 = betalisti[[2]]
  betafdPar2$lambda = lambda
  betalisti[[2]] = betafdPar2
  fRegi          = fRegress.CV(annualprec, templist, betalisti)
  SSE.CV[ilam]   = fRegi$SSE.CV
}
plot(loglam, SSE.CV, type="b", lwd=2,
     xlab="log smoothing parameter lambda",
     ylab="Cross-validation score", cex=1.2)

# Confidence Intervals
resid=annualprec-annualprechat
SigmaE.=sum(resid^2)/(35-fRegressList$df)
SigmaE=SigmaE.*diag(rep(1,35))
y2cMap=tempSmooth$y2cMap
stderrList=fRegress.stderr(fRegressList,y2cMap,SigmaE)

betafdPar = betaestlist2[[2]]
betafd = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd = betastderrList[[2]]
plot(betafd, xlab="Day",
     ylab="Temperature Reg. Coeff.",
     ylim=c(-6e-4,1.2e-03), lwd=2)
lines(betafd+2*betastderrfd, lty=2, lwd=1)
lines(betafd-2*betastderrfd, lty=2, lwd=1)

#Regression on FPCA
daybasis365   = create.fourier.basis(c(0, 365), 365)
lambda        = 1e6
tempfdPar365  = fdPar(daybasis365, harmaccelLfd, lambda)
tempSmooth365 = smooth.basis(day.5, daily$tempav,
                             tempfdPar365)
tempfd = tempSmooth365$fd

lambda    = 1e0
tempfdPar = fdPar(daybasis365, harmaccelLfd, lambda)
temppca   = pca.fd(tempfd, 4, tempfdPar)
harmonics = temppca$harmonics

pcamodel = lm(annualprec~temppca$scores)
pcacoefs = summary(pcamodel)$coef
betafd   = pcacoefs[2,1]*harmonics[1] + pcacoefs[3,1]*harmonics[2] +
  pcacoefs[4,1]*harmonics[3]
coefvar  = pcacoefs[,2]^2
betavar  = coefvar[2]*harmonics[1]^2 + coefvar[3]*harmonics[2]^2 +
  coefvar[4]*harmonics[3]^2

plot(betafd, xlab="Day", ylab="Regression Coef.",
     ylim=c(-6e-4,1.2e-03), lwd=2)
lines(betafd+2*sqrt(betavar), lty=2, lwd=1)
lines(betafd-2*sqrt(betavar), lty=2, lwd=1)


#CH10
#Setup Z
regions=unique(CanadianWeather$region)
p=length(regions)+1
regionList=vector('list',p)
regionList[[1]]=c(rep(1,35),0)
for(j in 2:p){
  xj=CanadianWeather$region==regions[j-1]
  regionList[[j]]=c(xj,1)
}
# fit tempfd
Lcoef        = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))
tempbasis    = create.fourier.basis(c(0, 365), 65)
lambda       = 1e6
tempfdPar65  = fdPar(tempbasis, harmaccelLfd, lambda)
dayOfYearShifted = c(182:365, 1:181)
tempShifted  = daily$tempav[dayOfYearShifted, ]
tempSmooth65 = smooth.basis(day.5, tempShifted, tempfdPar65)
tempfd       = tempSmooth65$fd
#add one more y of 0
coef=tempfd$coefs
coef36=cbind(coef,matrix(0,65,1))
temp36fd=fd(coef36,tempbasis,tempfd$fdnames)
#create functional parameter objects for each of the coefficient functions
betabasis=create.fourier.basis(c(0,365),11)
betafdPar=fdPar(betabasis)
betaList=vector('list',p)
for(j in 1:p)betaList[[j]]=betafdPar
fRegressList=fRegress(y=temp36fd,xfdlist=regionList,
                      betalist = betaList)
#Extract results
betaestList=fRegressList$betaestlist
regionFit=fRegressList$yhatfd
plot(temp36fd[1])
plot(regionFit)
regions=c('Canada',regions)
par(mfrow=c(2,3),cex=1)
for(j in 1:p) plot(betaestList[[j]]$fd,
                   lwd=2,xlab='Day',ylab='', 
                   main=regions[j])
plot(regionFit, lwd=2, col=1, lty=1,
     xlab="Day", ylab="",
     main="Prediction")

#Trends in Seabird Populations on Kodiak Island
# process data from seabird data.frame
sites = c('Uganik', 'Uyak')
sel   = seabird$Bay %in% sites
UU    = seabird[sel,]
# Drop 2 species with many NAs
NAs       = sapply(UU, function(x)sum(is.na(x)))
NAs.      = which(NAs > 2)
birdindex = (1:15)[-NAs.]
birds     = names(UU)[birdindex]
#  Compute mean counts taken over both sites and transects
meanCounts = matrix(NA, 20, 13)
dimnames(meanCounts) = list(1986:2005, birds)
for(i in 1:20){
  sel = (UU$Year == rownames(meanCounts)[i])
  meanCounts[i, ] = sapply(UU[sel, birds], mean, na.rm=TRUE)
}
selYear   = !is.na(meanCounts[,1])
logCounts = log10(meanCounts[selYear,])

#  time vectors in years and in indices in 1:20
yearObs  = as.numeric(rownames(logCounts))
yearCode = (1:20)[selYear]
shellfishindex = c(1,2,5,6,12,13)
fishindex      = (1:13)[-shellfishindex]

# Figure 10.2
ylim = range(logCounts)
op = par(mfrow=c(2,1), mar=c(2, 4, 4, 1)+.1)

matplot(yearObs, logCounts[, shellfishindex], xlab='', ylab='',
        ylim=ylim, main='Shellfish Diet', type='b', col=1)
meanShellfish = apply(meanCounts[, shellfishindex], 1, mean)
lines(yearObs, log10(meanShellfish[!is.na(meanShellfish)]), lwd=3)
abline(h=0, lty='dotted')

matplot(yearObs, logCounts[, fishindex], xlab='', ylab='',
        ylim=ylim, main='Fish Diet', type='b', col=1)
meanFish = apply(meanCounts[, shellfishindex], 1, mean)
lines(yearObs, log10(meanFish[!is.na(meanFish)]), lwd=3)
abline(h=0, lty='dotted')

par(op)

#  Compute mean counts taken over transects only within sites
#  so we have 2 observations for each bird species each year.
#  Two of these counts are zero, and are replaced by 1/(2*n)

meanCounts2 = matrix(NA, 20, 26)
for(i in 1:20) for (j in 1:2) {
  sel = (UU$Year == rownames(meanCounts)[i] & as.character(UU$Bay) == sites[j])
  meanCountsij = sapply(UU[sel, birds], mean, na.rm=TRUE)
  n = sum(sel)
  if (n > 0) {
    meanCountsij[meanCountsij == 0] = 1/(2*n)
  }
  meanCounts2[i,(j-1)*13+(1:13)] = meanCountsij
}

selYear2   = !is.na(meanCounts2[, 1])
yearCode  = (1:20)[selYear2]
all.equal(yearCode, c(1:12, 14:20))

logCounts2 = log10(meanCounts2[selYear2,])

#  Represent log mean counts exactly with a polygonal basis
birdbasis = create.polygonal.basis(yearCode)
birdlist2 = smooth.basis(yearCode, logCounts2, birdbasis)
birdfd2 = birdlist2$fd

#Define variable for diet effect
foodindex=c(1,2,5,6,12,13)
foodvarbl=(2*rep(1:13 %in% foodindex,2)-1)
#Define variable for bird species
birdvarbl=diag(rep(1,13))
birdvarbl=rbind(birdvarbl,birdvarbl)
Zmat0=matrix(0,26,15)
Zmat0[,1] = 1
Zmat0[,2] = foodvarbl
Zmat0[,3:15] = birdvarbl

#Remove singularity or add constraints
Zmat=rbind(Zmat0, matrix(0,2,15))
fishindex=(1:13)[-foodindex]
Zmat[27,foodindex+2]=1
Zmat[28,fishindex+2]=1

#  Two extra dummy observations are added to the functional data
#  object for log counts, and two additional rows are added to
#  the design matrix to force the bird effects within each diet
#  group to equal 0.

birdfd3 = birdfd2
birdfd3$coefs = cbind(birdfd3$coefs, matrix(0,19,2))

p = 15
xfdlist = vector("list",p)
names(xfdlist) = c("const", "diet", birds)
betalist = xfdlist
for (j in 1:p) xfdlist[[j]] = Zmat[,j]

#  set up the functional parameter object for (the regression fns.
#  use cubic b-spline basis for intercept and food coefficients
betabasis1 = create.bspline.basis(c(1,20),21,4,yearCode)
Lfdobj1    = int2Lfd(2);
Rmat1      = eval.penalty(betabasis1, Lfdobj1)
lambda1    = 10
betafdPar1 = fdPar(betabasis1,Lfdobj1,lambda1,TRUE,Rmat1)
betalist[[1]] = betafdPar1
betalist[[2]] = betafdPar1
betabasis2 = create.constant.basis(c(1,20))
betafdPar2 = fdPar(betabasis2)
for (j in 3:15) betalist[[j]] = betafdPar2

birdRegress = fRegress(birdfd3, xfdlist, betalist)
betaestlist = birdRegress$betaestlist

op = par(mfrow=c(3,1))
plot(betaestlist$const$fd)
plot(betaestlist$diet$fd)
plot(betaestlist[[3]]$fd)
par(op)

# Choosing smoothing parameter by CV
#  Choose the level of smoothing by minimizing cross-validated
#  error sums of squares.

loglam = seq(-2,4,0.25)
SSE.CV = rep(0,length(loglam))
betafdPari = betafdPar1
for(i in 1:length(loglam)){
  print(loglam[i])
  betafdPari$lambda = 10^loglam[i]
  betalisti = betalist
  for (j in 1:2) betalisti[[j]] = betafdPari
  CVi = fRegress.CV(birdfd3, xfdlist, betalisti, CVobs=1:26)
  SSE.CV[i] = CVi$SSE.CV
}
#fig 10.4
plot(loglam,SSE.CV,type='b',cex.lab=1.5,cex.axis=1.5,lwd=2,
     xlab='log smoothing parameter',ylab='cross validated sum of squares')


#Confidence interval
birdYhatfdobj=birdRegress$yhatfdobj
birdYhatmat = eval.fd(yearCode, birdYhatfdobj$fd[1:26])
rmatb   = logCounts2 - birdYhatmat
SigmaEb = var(t(rmatb))
y2cMap.bird = birdlist2$y2cMap
birdStderrList = fRegress.stderr(birdRegress, y2cMap.bird,
                                 SigmaEb)
birdBeta.sdList = birdStderrList$betastderrlist

op = par(mfrow=c(2,1))
plotbeta(betaestlist[1:3], birdBeta.sdList[1:3])

par(op)






