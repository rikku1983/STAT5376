##Test lambda
library(fda)
daytime=(1:365)-0.5
JJindex=c(182:365,1:181)
tempmat=daily$tempav[JJindex,]
tempbasis=create.fourier.basis(c(0,365), 365)
harmaccelLfd = vec2Lfd(c(0,(2*pi/365)^2,1), c(0, 365))
tempfdPar=fdPar(tempbasis, Lfdobj=harmaccelLfd, lambda=10e2)

tempfd=smooth.basis(daytime,tempmat,tempfdPar)$fd
tempfd$fdnames=list("Day (July 2 to June 30)",
                    "Weather Station",
                    "Mean temperature (deg. C)")
plot(tempfd, col=1,lty=1)#fig4.1

plot(daytime, tempmat[,1],'l')


loglam=seq(-8,8,2)
par(mfrow=c(3,3))
for(i in 1:length(loglam)){
  tempfdPari=fdPar(tempbasis, Lfdobj=harmaccelLfd, 
                  lambda=10^loglam[i])
  tempfdi=smooth.basis(daytime,tempmat[,1],tempfdPari)
  plot(tempfdi$fd,main=paste('lambda=', 10^loglam[i],', df=' ,tempfdi$df,''))
}

loglam=seq(-8,8,2)
par(mfrow=c(3,3))
for(i in 1:length(loglam)){
  tempfdPari=fdPar(tempbasis, Lfdobj = 2, 
                  lambda=10^loglam[i])
  tempfdi=smooth.basis(daytime,tempmat[,1],tempfdPari)
  plot(tempfdi$fd,main=paste('lambda=', 10^loglam[i],', df=' ,tempfdi$df,''))
}




tempbasis2=create.bspline.basis(c(0,365), nbasis=365, norder=5)
for(i in 1:length(loglam)){
  tempfdPari=fdPar(tempbasis2, Lfdobj = 2, 
                   lambda=10^loglam[i])
  tempfdi=smooth.basis(daytime,tempmat[,1],tempfdPari)
  plot(tempfdi$fd,main=paste('lambda=', 10^loglam[i],', df=' ,tempfdi$df,''))
}

#fig 5.1 GCV choosing lambda
ageRng  = c(1,18)
age     = growth$age
heightmat  = growth$hgtf
norder = 6
nbasis = length(age) + norder - 2
heightbasis = create.bspline.basis(ageRng, nbasis, norder, age)

loglam = seq(-6, 0, 0.25)
Gcvsave = rep(NA, length(loglam))
names(Gcvsave) = loglam
Dfsave = Gcvsave
for(i in 1:length(loglam)){
  hgtfdPari = fdPar(heightbasis, Lfdobj=4, 10^loglam[i])
  hgtSm.i = smooth.basis(age, heightmat, hgtfdPari)
  Gcvsave[i] = sum(hgtSm.i$gcv)
  Dfsave[i] = hgtSm.i$df
}

# Figure 5.1.

plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(GCV(lambda)), lwd=2 )


# Smoothing case study
yearRng   = c(0,365)
dayOfYearShifted = c(182:365, 1:181)
logprecav = CanadianWeather$dailyAv[dayOfYearShifted, , 'log10precip']
nbasis   = 365
daybasis = create.fourier.basis(yearRng, nbasis)
Lcoef        = c(0,(2*pi/diff(yearRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, yearRng)
loglam        = seq(4,9,0.25)
nlam          = length(loglam)
dfsave        = rep(NA,nlam)
names(dfsave) = loglam
gcvsave       = dfsave
for (ilam in 1:nlam) {
  cat(paste('log10 lambda =',loglam[ilam],'\n'))
  lambda        = 10^loglam[ilam]
  fdParobj      = fdPar(daybasis, harmaccelLfd, lambda)
  smoothlist    = smooth.basis(day.5, logprecav,
                               fdParobj)
  dfsave[ilam]  = smoothlist$df
  gcvsave[ilam] = sum(smoothlist$gcv)
}

# Figure 5.2.
plot(day.5,logprecav[,1],'l')
plot(loglam, gcvsave, type='b', lwd=2, ylab='GCV Criterion',
     xlab=expression(log[10](lambda)) )

lambda = 1e6
fdParobj = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fit = smooth.basis(day.5, logprecav, fdParobj)
logprec.fd = logprec.fit$fd
fdnames = list("Day (July 1 to June 30)",
                   "Weather Station" = CanadianWeather$place,
                   "Log 10 Precipitation (mm)")
logprec.fd$fdnames = fdnames

#  plot the functional data object
plot(logprec.fd)
plotfit.fd(logprecav[,1], day.5, logprec.fd[1], lwd=1)

logprecmat=eval.fd(day.5, logprec.fd)
logprecres=logprecav-logprecmat
logprecvar1 = apply(logprecres^2, 1, sum)/35
logprecvar2=apply(logprecres^2,2,sum)/(365-12)
plot(sqrt(logprecvar2), xlab='Station Number',
     ylab='Standard Deviation across Day')
rt = which(CanadianWeather$place %in%
             c("Winnipeg", 'Regina', 'Churchill', 'Montreal', 'St. Johns'))
lft = which(CanadianWeather$place %in%
              c('Yellowknife', 'Resolute', 'Vancouver', 'Iqaluit',
                'Pr. George', 'Pr. Rupert') )
below = which(CanadianWeather$place %in% 'Edmonton')
top = which(CanadianWeather$place %in% 'Halifax')

text(rt, sqrt(logprecvar2[rt]), labels=CanadianWeather$place[rt],
     pos=4)
text(lft, sqrt(logprecvar2[lft]), labels=CanadianWeather$place[lft],
     pos=2)
text(below, sqrt(logprecvar2[below]), labels=CanadianWeather$place[below],
     pos=1)
text(top, sqrt(logprecvar2[top]), labels=CanadianWeather$place[top],
     pos=3)


logstddev.fd = smooth.basis(day.5,log(logprecvar1)/2, fdParobj)$fd
logprecvar1fit = exp(eval.fd(day.5, logstddev.fd))
plot(day.5, sqrt(logprecvar1), xlab='Day',ylab='Standard seviation across stations')
lines(day.5, logprecvar1fit, lwd=2)

#Bivariate covariance function
logprecvar.bifd=var.fd(logprec.fd)
weektime=seq(0,365,length=53)
logprecvar_mat=eval.bifd(weektime, weektime, logprecvar.bifd)
#persp(weektime,weektime,logprecvar_mat, theta=-45, phi=25,r=3,expand=0.5,)

library(reshape2)
df=melt(logprecvar_mat)
write.csv(df,'prec_var_surface.csv2')

# 
# names(df)=c('weektime1','weektime2','Covariance')
# library(ggplot2)
# p=ggplot(df, aes(x=weektime1, y=weektime2))+geom_tile(aes(fill=Covariance))
# library(plotly)
# Sys.setenv("plotly_username" = "rikku1983")
# Sys.setenv("plotly_api_key" = "ph0KqUyMlpTvBwYxfzqD")
# 
# plotly_POST(p, filename = "Covar_surface")
# ggplotly(p)
# # 
# # dsamp <- diamonds[sample(nrow(diamonds), 1000), ]
# # qplot(carat, price, data=dsamp, colour=clarity)
# # ggplotly()
# # set.seed(100)
# # d <- diamonds[sample(nrow(diamonds), 1000), ]
# # 
# # p <- ggplot(data = d, aes(x = carat, y = price)) +
# #   geom_point(aes(text = paste("Clarity:", clarity)), size = 4) +
# #   geom_smooth(aes(colour = cut, fill = cut)) + facet_wrap(~ cut)
# # 
# # (gg <- ggplotly(p))
# # 
# # py <- plotly(username="rikku1983", key="8SBOMDYBiI71EF7tycTr")
# # ggplotly(p)


dayvec = seq(0,365,len=101)
xivec = exp(20*cos(2*pi*(dayvec-197)/365))
xibasis = create.bspline.basis(c(0,365),13)
xifd = smooth.basis(dayvec, xivec, xibasis)$fd
tempbasis = create.fourier.basis(c(0,365),65)
precbasis = create.fourier.basis(c(0,365),365)
tempLmat = inprod(tempbasis, xifd)
precLmat = inprod(precbasis, xifd)
test=inprod(logprec.fd,xifd)
test
test2=inprod(tempfd,xifd)
plot(test2)

#Calculate CL
#  smooth data with lambda that minimizes GCV getting
#  all of the output up to matrix y2cMap
yearRng      = c(0,365)
daybasis     = create.fourier.basis(yearRng, 365)
Lcoef        = c(0,(2*pi/diff(yearRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, yearRng)
lambda     = 1e6
fdParobj   = fdPar(daybasis, harmaccelLfd, lambda)

logprecList = smooth.basis(day.5, logprecav, fdParobj)
logprec.fd  = logprecList$fd
fdnames    = list("Day (July 1 to June 30)",
                  "Weather Station" = CanadianWeather$place,
                  "Log10 Precipitation (mm)")
logprec.fd$fdnames = fdnames

#  compute the residual matrix and variance vector
logprecmat  = eval.fd(day.5, logprec.fd)
logprecres  = logprecav - logprecmat
logprecvar  = apply(logprecres^2, 1, sum)/(35-1)

#  smooth log variance vector
lambda      = 1e8
resfdParobj = fdPar(daybasis, harmaccelLfd, lambda)
logvar.fd   = smooth.basis(day.5, log(logprecvar), resfdParobj)$fd

#  evaluate the exponentiated log variance vector and
#  set up diagonal error variance matrix SigmaE
varvec      = exp(eval.fd(day.5, logvar.fd))
SigmaE      = diag(as.vector(varvec))

#  compute variance covariance matrix for fit
y2cMap        = logprecList$y2cMap
c2rMap        = eval.basis(day.5, daybasis)
Sigmayhat     = c2rMap %*% y2cMap %*% SigmaE %*%
  t(y2cMap) %*% t(c2rMap)

#  extract standard error function for yhat

logprec.stderr= sqrt(diag(Sigmayhat))

#  plot Figure 6.6

logprec29 = eval.fd(day.5, logprec.fd[29])

plot(logprec.fd[29], lwd=2, ylim=c(0.2, 1.3))
lines(day.5, logprec29 + 2*logprec.stderr,
      lty=2, lwd=2)
lines(day.5, logprec29 - 2*logprec.stderr,
      lty=2, lwd=2)
points(day.5, logprecav[,29])


#PCA
tempfdpca=pca.fd(tempfd,3, centerfns = F)
s=tempfdpca$scores
basis=eval.fd(day.5, tempfdpca$harmonics)
sim=basis %*% t(s)
matplot(sim)
plot(tempfd)

#CH8 Registration

#CH9 
annualprec = log10(apply(daily$precav,2,sum))
tempbasis65=create.fourier.basis(c(0,365), 65)
tempSmooth65 = smooth.basis(day.5, daily$tempav, tempbasis65)
tempfd65     = tempSmooth65$fd

templist = vector("list",2)
templist[[1]] = rep(1,35)
templist[[2]] = tempfd65

conbasis = create.constant.basis(c(0,365))
betabasis = create.fourier.basis(c(0,365),5)
betalist = vector("list",2)
betalist[[1]] = conbasis
betalist[[2]] = betabasis

fRegressList = fRegress(annualprec, templist, betalist)
betaestlist = fRegressList$betaestlist
tempbetafd = betaestlist[[2]]$fd
plot(tempbetafd, xlab="Day",ylab="Beta for temperature")

yhat=fRegressList$yhatfdobj
plot(1:35,annualprec,'l',col='black', lwd=2,ylab='Annual precipitation', xlab='Stations')
lines(1:35, yhat, lwd=2,col='red')
legend(2,2.4, c('Real', 'fitted'),lty=c(1,1), lwd=c(2.5,2.5),col=c('black','red'))

annualprechat1 = fRegressList$yhatfdobj
annualprecres1 = annualprec - annualprechat1
SSE1.1 = sum(annualprecres1^2)
SSE0 = sum((annualprec - mean(annualprec))^2)
SSR=sum((annualprechat1-mean(annualprec))^2)
RSQ1 = (SSE0-SSE1.1)/SSE0
Fratio1 = ((SSE0-SSE1.1)/5)/(SSE1.1/29)

# regression using roughness penalty
Lcoef = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))

# refit with 35 terms rather than 5 in the fourier basis

betabasis35 = create.fourier.basis(c(0, 365), 35)
lambda = 10^12.5
betafdPar. = fdPar(betabasis35, harmaccelLfd, lambda)
betalist2 = betalist1
betalist2[[2]] = betafdPar.


annPrecTemp = fRegress(annualprec, templist, betalist2)
betaestlist2 = annPrecTemp$betaestlist
plot(betaestlist2[[2]]$fd, xlab='Day', ylab='Beta for temperature')
betaestlist2[[1]]$fd
annPrecTemp$df
yhat2=annPrecTemp$yhatfdobj
plot(1:35,annualprec,'l',col='black', lwd=2,ylab='Annual precipitation', xlab='Stations')
lines(1:35, yhat2, lwd=2,col='red')
legend(2,2.4, c('Real', 'fitted'),lty=c(1,1), lwd=c(2.5,2.5),col=c('black','red'))


annualprechat2 = annPrecTemp$yhatfdobj
print(annPrecTemp$df)
SSE1.2 = sum((annualprec-annualprechat2)^2)
(RSQ2 = (SSE0 - SSE1.2)/SSE0)
(Fratio2 = ((SSE0-SSE1.2)/3.7)/(SSE1.2/30.3))

#Residual plot
plot(annualprechat2, annualprec, lwd=2, col='red', xlab='Predectied annual precipitation')
abline(lm(annualprec~annualprechat2), lty='dashed', lwd=2)
points(yhat,annualprec,lwd=2, col='blue')
legend(2.25,3.4,c('Model1','Model2'), pch=c(1,1), lwd=c(2,2),col=c('blue','red'))


#Confidence interval
#Model 1
resid   = annualprec - yhat
SigmaE. = sum(resid^2)/(35-fRegressList$df)
SigmaE  = SigmaE.*diag(rep(1,35))
#y2cMap  = tempSmooth65$y2cMap
stderrList = fRegress.stderr(fRegressList, diag(rep(1,35)), SigmaE)

betafdPar      = betaestlist[[2]]
betafd         = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd   = betastderrList[[2]]

plot(betafd, xlab="Day", ylab="Temperature Reg. Coeff.",
     ylim=c(-2.3e-3,2.2e-03), lwd=2)
lines(betafd+2*betastderrfd, lty=2, lwd=1)
lines(betafd-2*betastderrfd, lty=2, lwd=1)

#Model 2
resid   = annualprec - annualprechat2
SigmaE. = sum(resid^2)/(35-annPrecTemp$df)
SigmaE  = SigmaE.*diag(rep(1,35))
#y2cMap  = tempSmooth65$y2cMap
stderrList = fRegress.stderr(annPrecTemp, diag(rep(1,35)), SigmaE)

betafdPar      = betaestlist2[[2]]
betafd         = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd   = betastderrList[[2]]

plot(betafd, xlab="Day", ylab="Temperature Reg. Coeff.",
     ylim=c(-6e-4,1.2e-03), lwd=2)
lines(betafd+2*betastderrfd, lty=2, lwd=1)
lines(betafd-2*betastderrfd, lty=2, lwd=1)

F.res = Fperm.fd(annualprec, templist, betalist2)
F.res$Fobs
F.res$qval


# CH10

regions. = unique(CanadianWeather$region)
p = length(regions.) + 1
regionList = vector("list", p)
names(regionList) = c('Canada', regions.)
regionList[[1]] = c(rep(1,35),0)
for (j in 2:p) {
  xj = (CanadianWeather$region == regions.[j-1])
  regionList[[j]]= c(xj,1)
}

w<- as.data.frame(CanadianWeather$coordinates)
w$region=CanadianWeather$region
library(ggmap)
map=get_map('Canada', zoom = 3)
ggmap(map) + geom_point(data=w, aes(x=-W.longitude, y=N.latitude, col=region), size=3)

Lcoef        = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))
tempbasis    = create.fourier.basis(c(0, 365), 65)
lambda       = 1e6
tempfdPar65  = fdPar(tempbasis, harmaccelLfd, lambda)
dayOfYearShifted = c(182:365, 1:181)
tempShifted  = daily$tempav[dayOfYearShifted, ]
tempSmooth65 = smooth.basis(day.5, tempShifted, tempfdPar65)
tempfd       = tempSmooth65$fd
coef = tempfd$coef
coef36 = cbind(coef,matrix(0,65,1))
temp36fd = fd(coef36,tempbasis,tempfd$fdnames)

betabasis = create.fourier.basis(c(0, 365), 11)
betafdPar = fdPar(betabasis)
betaList = vector("list",p)
names(betaList) = regions.
for (j in 1:p) betaList[[j]] = betafdPar

fRegressList= fRegress(temp36fd, regionList, betaList)
betaestList = fRegressList$betaestlist
regionFit = fRegressList$yhatfd
regions = c("Canada", regions.)

op = par(mfrow=c(2,3),cex=1)
for (j in 1:p) plot(betaestList[[j]]$fd, lwd=2,
                    xlab="Day (July 1 to June 30)",
                    ylab="", main=regions[j])
plot(regionFit$fd, lwd=2, col=1, lty=1,
     xlab="Day (July 1 to June 30)", ylab="", main="Prediction")
par(op)

fd=vector('list',4)
for (i in 1:4){
  fd[[i]]=mean(tempfd[CanadianWeather$region==regions.[i]])
}
par(mfrow=c(1,1))
plot(fd[[1]],lwd=2, ylim=c(-32,20))
plot(fd[[2]],lwd=2, add=T)
plot(fd[[3]],lwd=2, add=T)
plot(fd[[4]],lwd=2, add=T)

y2cmap2=tempSmooth65$y2cMap

#variance
yhatfdobj2=fRegressList$yhatfdobj$fd
yhatmat = eval.fd(day.5, yhatfdobj2)[,1:35]
ymat = tempShifted
rmatb = ymat - yhatmat
SigmaEb = var(t(rmatb))
y2cMap = tempSmooth65$y2cMap
stderrList = fRegress.stderr(fRegressList, y2cMap,SigmaEb)
str(stderrList)
betalist3=stderrList$betastderrlist
class(betalist3[[1]])
plotbeta(betaestList,betalist3)
