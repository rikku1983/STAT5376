#Warping function, taking 2 fd list objects, rfd should be defined on [0,1], ffdlist should be defined from 0
warp=function(fd, rfd){
  n=length(fd$fdnames$time)
  basisf=fd$basis
  frange=fd$basis$rangeval
  xf=seq(frange[1],frange[2],len=n)
  xfs=(xf-frange[1])/frange[2]
  newx=eval.fd(xfs,rfd)*frange[2]+frange[1]
  newx[1]=frange[1]
  newx[n]=frange[2]
  y=eval.fd(newx, fd)
  newffd=smooth.basis(xf,y,basisf)
  return(newffd)
}

#Function looking for neighbors
neighbors=function(nn){
  neighbs=matrix(1,1,2)
  slope=1
  for(i in 1:nn){
    for(j in 1:nn){
      if(!((i/j) %in% slope)){
        neighbs=rbind(neighbs,c(i,j))
        slope=c(slope,i/j)
      }
    }
  }
  return(neighbs)
}

#Random warping function object generator, xgrid is a vector from 0. nn is scope of neighbors with default value
rndwarper=function(xgrid, nn=ceiling(length(xgrid)/20)){
  #define possible next points with different slope
  n=length(xgrid)-1
  nei=neighbors(nn)
  m=min(nn,n-1)
  slopelim=c(1/(m),(m))
  curpnt=c(0,0)
  path=curpnt
  while(curpnt[1] < n){
    #find possible point to go
    posspnt=t(t(nei)+curpnt)
    impo=vector()
    for(i in 1:nrow(posspnt)){
      tstslp=(n-posspnt[i,1])/(n-posspnt[i,2])
      if(is.na(tstslp)){tstslp=999999}
      if(posspnt[i,1]>n | posspnt[i,2]>n | tstslp < slopelim[1] | tstslp>slopelim[2]){
        impo[i]=FALSE
        if(posspnt[i,1]==n & posspnt[i,2]==n){impo[i]=TRUE}
      }
      else{impo[i]=TRUE}
    }
    posspnt=posspnt[impo,]
    if(length(posspnt)==2){
      curpnt=posspnt
    }
    else{
      idx=sample(1:(length(posspnt)/2),1)
      curpnt=posspnt[idx,]
    }
    path=rbind(path,curpnt)
  }
  x=path[,1]/max(xgrid)
  y=path[,2]/max(xgrid)
  #interpolate
  rbasis=create.bspline.basis(c(0,1),norder = 2, breaks = x)
  rfd=Data2fd(x, y, rbasis)
  return(rfd)
}


# Q function
f2q=function(fd,x){
  basisf = fd$basis
  frange = fd$basis$rangeval
  d1y = eval.fd(x, fd, 1)
  qv=sign(d1y)*sqrt(abs(d1y))
  q=smooth.basis(x,qv,basisf)$fd
}

# Path to warper function
pth2r=function(path){
  path=path[order(path[,1]),]
  # path=path[-1,]
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

# Dynamic programming function to find path with lowest cost, taking discrete valued functions
sldp2=function(f,g){
  c=Inf
  fn=length(f)
  n=101
  m=length(g)
  xx1=(0:fn-1)/(fn-1)
  xx2=(0:m-1)/(m-1)
  E=matrix(0,n,n)
  E[1,]=c
  E[,1]=c
  E[1,1]=0
  path=array(dim=c(n,n,2))
  v=neighbors(6)
  for(i in 2:n){
    for(j in 2:n){
      # print(paste('i=',i))
      # print(paste('j=',j))
      cande=as.numeric()
      for(r in 1:nrow(v)){
        # print(paste('r=',r))
        k=i-v[r,1]
        l=j-v[r,2]
        if(k>0 & l>0){
          cande[r]=E[k,l]+energySRSF2(f,g,k,l,i,j)
          # print(energySRSF2(f,g,k,l,i,j))
          }
          else{
            cande[r]=c
          }
        }
        # print(cande)
        minidx=which.min(cande)
        E[i,j]=cande[minidx]
        path[i,j,1]=i-v[minidx,1]
        path[i,j,2]=j-v[minidx,2]
      }
    }
    pth=matrix(c(n,n),1,2)
    nextpnt=c(n,n)
    while(nextpnt[1]!=1 & nextpnt[2]!=1){
      curpnt=pth[nrow(pth),]
      nextpnt=path[curpnt[1],curpnt[2],]
      pth=rbind(pth, nextpnt)
    }
    return(list(path, pth, E))
}


# Energy function taking discrete functions  
energySRSF2=function(f,g,k,l,i,j){
    fn=length(f)
    m=length(g)
    slope=(j-l)/(i-k)
    fidx=round(k*fn/101):round(i*fn/101)
    gidx=round(l*m/101):round(j*m/101)
    cost=sum((f[fidx]-g[gidx]*sqrt(slope))^2)/n
    return(cost)
}

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