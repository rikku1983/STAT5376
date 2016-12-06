library(fda)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
hw=handwrit

#Simulate function
par(mfrow=c(1,2))
x=seq(0, 10,len=1001)
y1=dnorm(x, mean=3, sd=0.3)+dnorm(x, mean=4, sd=0.3)
y2=dnorm(x, mean=5, sd=0.3)+dnorm(x, mean=6, sd=0.3)
bas=create.bspline.basis(range(x), norder=5, nbasis=50)
fd1=smooth.basis(x,y1,bas)$fd
fd2=smooth.basis(x,y2,bas)$fd

q1=f2q(fd1,x)
q2=f2q(fd2,x)

result1 = sldp(fd1,fd2)
pth = result[[2]]
# plot(pth[,1],pth[,2])
r1=pth2r(pth)
# plot(r)
nfd2=warp(fd2,r)

result2 = sldp(q1,q2)
pth = result[[2]]
# plot(pth[,1],pth[,2])
r2=pth2r(pth)
# plot(r)
nfd2q=warp(fd2,r2)

plot(fd1, ylim=c(0,1.5))
plot(fd2, col='red', add=T)
plot(nfd2, col='blue', add=T)

plot(fd1,ylim=c(0,1.5))
plot(fd2, col='red', add=T)
plot(nfd2q, col='blue', add=T)


plot(q1)

r1=pth2r(result[[2]])
gfdw=warp(d1list[[2]], r1, 601)

# Path to warper function
pth2r=function(path){
  path=path[order(path[,1]),]
  path=path[-1,]
  path=path-1
  x=path[,1]/max(path[,1])
  y=path[,2]/max(path[,2])
  #interpolate
  rbasis=create.bspline.basis(c(0,1),norder = 2, breaks = x)
  rfd=Data2fd(x, y, rbasis)
}



# My plot function

fdaplot=function(x,y,colr){
  df=data.frame(x=x,y=y,colr=colr)
  p=ggplot(aes(x,y,colour=colr), data=df)+geom_point()+scale_colour_gradientn(colors=c('black','blue','green','yellow','red'))
  return(p)
}


tx=hw[,1,1]
xbas=create.bspline.basis(c(0,1400), norder = 6, nbasis = 100)
txfd=smooth.basis(0:1400, tx, xbas)
plot(txfd)
lines(0:1400, tx, col='red')
#100 basis is enough
rfd=rndwarper(0:100, 2)
newtxfdlist=warp(txfd,rfd)
plot(newtxfdlist$fd, add=TRUE)

par(mfrow=c(1,2))
plot(hw[,1,1],hw[,1,2])
plot(eval.fd(0:1400, newtxfdlist$fd), hw[,1,2])

tyfd=smooth.basis(0:1400, hw[,1,2], xbas)
newtyfdlist=warp(tyfd,rfd)
plot(tyfd$fd)
plot(newtyfdlist$fd, col='red',add=T)

plot(eval.fd(0:1400, newtxfdlist$fd), eval.fd(0:1400, newtyfdlist$fd))
mypal=colorRampPalette(c('black','blue','green','yellow','red'))(1401)

fdaplot(hw[,1,1], hw[,1,2], 0:1400)
fdaplot(eval.fd(0:1400, newtxfdlist$fd), eval.fd(0:1400, newtyfdlist$fd), 0:1400)

plot(x=hw[,1,1],y=hw[,1,2], 'l', col=mycol[1])
for(i in 2:dim(hw)[2]){
  lines(x=hw[,i,1],y=hw[,i,2],col=mycol[i])
}
plot(x=1:dim(hw)[1],y=hw[,1,1], 'l', col=mycol[1])
for(i in 2:dim(hw)[2]){
  lines(x=1:dim(hw)[1],y=hw[,i,1],col=mycol[i])
}

plot(x=1:dim(hw)[1],y=hw[,1,2], 'l', col=mycol[1])
for(i in 2:dim(hw)[2]){
  lines(x=1:dim(hw)[1],y=hw[,i,2],col=mycol[i])
}
p=ggplot(aes(x,y,colour=col), data=df)+geom_point()+scale_colour_gradientn(colors=c('black','blue','green','yellow','red'))


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
plot(d1)

fdaplot(ssc[,1,1], ssc[,1,2], 0:600)

# Function to return writing speed
wspeed=function(t,x,y,nbas){
  bas=create.bspline.basis(range(t), norder=6, nbasis=nbas)
  xfd=smooth.basis(t, x, bas)
  yfd=smooth.basis(t, y, bas)
  xfd1=deriv.fd(xfd$fd, 1)
  yfd1=deriv.fd(yfd$fd, 1)
  d1=sqrt((xfd1+1)^2+(yfd1+1)^2)
  return(d1)
}
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
# hwd1list=vector('list',20)
# for(i in 1:20){
#   hwd1list[[i]]=wspeed(0:1400, hw[,i,1], hw[,i,2], 100)
#   if(i==1){
#     plot(hwd1list[[i]])
#     next
#   }
#   plot(hwd1list[[i]],add=T)
# }

# Convert all speed function to q
qlist=vector('list', 50)
for(i in 1:50){
  qlist[[i]] = f2q(d1list[[i]], 0:600)
}
saveRDS(qlist,'qlist.rds')


re=sldp(qlist[[1]]$fd, qlist[[2]]$fd)

# Energy function
energySRSF = function(fd, gd, k, l, i, j, n){
  slope=(j-l)/(i-k)
  fx=(seq(k,i,len=100)-1)*fd$basis$rangeval[2]/n
  gx=(seq(l,j,len=100)-1)*gd$basis$rangeval[2]/n
  fy=eval.fd(fx,fd)
  gy=eval.fd(gx,gd)
  cost=sum((fy-gy*sqrt(slope))^2)
  return(cost)
}

# DP function
sldp = function(fd, gd){
  c = Inf
  n = 10
  E = matrix(0,n+1,n+1)
  E[1,] = c
  E[,1] = c
  E[1,1] = 0
  path=array(dim=c(n+1,n+1,2))
  path[1,1,] = c(0,0)
  v = neighbors(10)
  for(i in 2:(n+1)){
    for(j in 2:(n+1)){
      candE = as.numeric()
      for(r in 1:nrow(v)){
        k = i-v[r,1]
        l = j-v[r,2]
        if(k>0 & l>0){
          candE[r] = E[k,l] + energySRSF(fd,gd,k,l,i,j,n)}
        else{
          candE[r] = c
        }
      }
      minidx = which.min(candE)
      E[i,j]=candE[minidx]
      path[i,j,1] = i-v[minidx,1]
      path[i,j,2] = j-v[minidx,2]
    }
  }
  pth=matrix(c(n+1,n+1),1,2)
  nextpnt=c(n+1,n+1)
  while(nextpnt[1]!=0 & nextpnt[2]!=0){
    curpnt=pth[nrow(pth),]
    nextpnt=path[curpnt[1],curpnt[2],]
    pth=rbind(pth, nextpnt)
  }
  return(list(path, pth, E))
}


### Test dynamic function for registration
n=1000
m=1000
xf=(0:n)/n
xg=(0:m)/m
f=dnorm(xf, mean=0.3, sd=0.03)#+dnorm(xf, mean=0.4, sd=0.04)
g=dnorm(xf, mean=0.5, sd=0.05)#+dnorm(xg, mean=0.6, sd=0.03)

re=sldp2(f,g)
r=pth2r(re[[2]])
plot(r)

bas=create.bspline.basis(range(xf), norder=5, nbasis=50)
ffd=smooth.basis(xf,f,bas)$fd
gfd=smooth.basis(xg,g,bas)$fd
ngfd=warp(gfd,r)

plot(ffd)
plot(gfd, col='red', add=T)
plot(ngfd, col='blue', add=T)
