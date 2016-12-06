library(gridExtra)
library(animation)
library(fda)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

# Load data
ssc=StatSciChinese
xbas=create.bspline.basis(c(0,600), norder = 6, nbasis = 70)
xfd=smooth.basis(0:600, ssc[,1,1], xbas)
plot(0:600, ssc[,1,1],'l', ylab='x coordinate', xlab='time')
plot(xfd$fd, col='red', add=T)

ybas=create.bspline.basis(c(0,600), norder = 6, nbasis = 70)
yfd=smooth.basis(0:600, ssc[,1,2], ybas)
plot(0:600, ssc[,1,2],'l', ylab='y coordinate', xlab='time')
plot(yfd$fd, col='red', add=T)

xfd1=deriv.fd(xfd$fd, 1)
plot(xfd1)
yfd1=deriv.fd(yfd$fd,1)
plot(yfd1)

d1=sqrt(xfd1^2+yfd1^2)
plot(d1, ylab='Speed of pen moving', main = 'Pen moving speed curve')

fdaplot(ssc[,1,1], ssc[,1,2], 0:600)

mycol=colorRampPalette(c('black','blue','green','yellow','red'))(50)
d1list=vector('list',50)
for(i in 1:50){
  d1list[[i]]=wspeed(0:600, ssc[,i,1], ssc[,i,2], 80)
  if(i==1){
    plot(d1list[[i]], col=mycol[i], ylab='adjusted pen speed')
    next
  }
  plot(d1list[[i]],add=T, col=mycol[i])
}
saveRDS(d1list,'d1list.rds')

# Convert all speed function to q
qlist=vector('list', 50)
for(i in 1:50){
  qlist[[i]] = f2q(d1list[[i]], 0:600)
}
saveRDS(qlist,'qlist.rds')

plot(qlist[[2]], ylab='q', main='SRSF representation')
mycol=colorRampPalette(c('black','blue','green','yellow','red'))(50)

# Calculate mean function
sumf = qlist[[1]]
plot(sumf, col=mycol[1])
for(i in 2:50){
  plot(qlist[[i]],col=mycol[i], add=T)
  sumf = sumf + qlist[[i]]
}
meanq = sumf*0.02
plot(meanq, lty=2, lwd=3, add=T)

# Registration
x=seq(0, 600, len=1001)
q1v=eval.fd(x,qlist[[1]])
dp=0
da=0
d1=d1list[[1]]
newd=vector('list',50)
newd[[1]]=d1
sumd=d1
plot(d1, col=mycol[1], ylab='Pen moving speed', main='Registered observations to 1st observation')
for(i in 2:50){
  q2v=eval.fd(x,qlist[[i]])
  # plot(x,q1v,'l')
  # lines(x,q2v,col='red')
  re=sldp2(q1v,q2v)
  r12=pth2r(re[[2]])
  d2=d1list[[i]]
  newd2=warp(d2,r12)$fd
  newd[[i]]=newd2
  sumd = sumd + newd2
  plot(newd2, col=mycol[i], add=T)
  # plot(d2, col='blue', add=T)
  # legend(x=250, y=1.425, c('f1(t)', 'f2(r(t))', 'f2(t)'), col=c('black', 'red', 'blue'), lty=1)
  
  #Calculate phase difference dp
  dp[i] = acos(sum(sqrt(eval.fd(seq(0,1,len=100)[2:99],r12,1)))/98)
  #Calculate amplitude difference da
  da[i] = sqrt(sum((eval.fd(seq(0,600,len=1000)[2:999],d1-newd2))^2)/998)
}
saveRDS(newd, 'newd.rds')
meand = sumd*0.02
saveRDS(meand, 'meand.rds')
plot(meand, lty=2, lwd=3, add=T)

re=sldp2(q1v,q1v)
r=pth2r(re[[2]])
plot(r)
acos(sum(sqrt(eval.fd(seq(0,1,len=100)[2:99],r,1)))/98)

aniplot=function(t,nth){
  df=data.frame(time=0:600, x=ssc[,nth,1], y=ssc[,nth,2], spd=eval.fd(0:600, meand))
  p1=ggplot(data=df[df$time<=t,], aes(x,y,colour=time))+geom_point()+ scale_x_continuous(limits=c(-0.18, 0.205)) + scale_y_continuous(limits=c(-0.1, 0.1))+scale_colour_gradientn(colors=c('black','blue','green','yellow','red'))
  p2=ggplot(data=df[df$time<=t,], aes(x=time, y=spd, colour=time))+geom_line(size=2)+scale_x_continuous(limits=c(0,600)) + scale_y_continuous(limits=c(1.405, 1.425))+scale_colour_gradientn(colors=c('black','blue','green','yellow','red'))
  grid.arrange(p1, p2, nrow=2)
}

trace.animate <- function(nth) {
  lapply(seq(0,600,length.out = 200), function(i) {
    aniplot(i,nth)
  })
}

#save all iterations into one GIF
saveGIF(trace.animate(1), interval = .2, movie.name="trace.gif")

install.packages("devtools")
library(devtools)

saveGIF({
  for (i in 1:10) plot(runif(10), ylim = 0:1)
})
