find_stack <- function (ijn,ndon,donor){
  for (i in 1:(ndon[ijn])){
    ijk<-donor[i,ijn]
    nstack<<-nstack+1
    stack[nstack]<<-ijk
    if (ndon[ijk]!=0) find_stack(ijk,ndon,donor)
  }
}

nx<-51
ny<-51
nn<-nx*ny
h<-array(0,nn)
rec<-array(0,nn)
le<-array(0,nn)
ndon<-array(0,nn)
donor<-array(0,c(8,nn))
stack<-array(0,nn)
a<-array(0,nn)

xl<-100000.
yl<-100000.
dx<-xl/(nx-1)
dy<-yl/(ny-1)

precipitation<-1.
kf<-1.e-4
m<-0.4
dt<-10000.
nstep<-100
u<-2.e-3
kd<-0

h<-runif(nn)

for (istep in 1:nstep){

for (ij in 1:nn){rec[ij]<-ij}

for (j in 2:(ny-1)){
  for (i in 2:(nx-1)){
    slopemax<-0.0000001
    ij<-i+(j-1)*nx
    for (jj in -1:1){
      for (ii in -1:1){
        iii<-i+ii
        jjj<-j+jj
        ijk<-iii+(jjj-1)*nx
        if (ijk!=ij){
          leng<-sqrt((dx*ii)^2+(dy*jj)^2)
          slope<-(h[ij]-h[ijk])/leng
          if (slope>slopemax){
            slopemax<-slope
            rec[ij]<-ijk
            le[ij]<-leng
          }
        }
      }
    }
  }
}

for (ij in 1:nn){ndon[ij]<-0}

for (ij in 1:nn){
  if (ij!=rec[ij]){
    ijr<-rec[ij]
    ndon[ijr]<-ndon[ijr]+1
    donor[ndon[ijr],ijr]<-ij
  }
}

nstack<-0
for (ij in 1:nn){
  if (ij==rec[ij]){
    nstack<-nstack+1
    stack[nstack]<-ij
    if (ndon[ij]!=0) find_stack(ij,ndon,donor)
  }
}

for (ij in 1:nn){a[ij]<-precipitation*dx*dy}

for (ij in rev(stack)) {
  if (rec[ij]!=ij){
    ijr<-rec[ij]
    a[ijr]<-a[ijr]+a[ij]
  }
}

for (j in 2:(ny-1)){
  for (i in 2:(nx-1)){
    ij<-i+(j-1)*nx
    h[ij]=h[ij]+u*dt
    }
  }

for (ij in stack){
  if (rec[ij]!=ij){
    ijr<-rec[ij]
    fact<-dt*kf*a[ij]^m/le[ij]
    h[ij]<-(h[ij]+fact*h[ijr])/(1.+fact)
  }
}

factx=dt*kd/dx^2
facty=dt*kd/dy^2
for (j in 2:(ny-1)){
  for (i in 2:(nx-1)){
    ij=i+(j-1)*nx
    imj=i-1+(j-1)*nx
    ipj=i+1+(j-1)*nx
    ijm=i+(j-2)*nx
    ijp=i+j*nx
    h[ij]=h[ij]+factx*(h[ipj]-2*h[ij]+h[imj])+facty*(h[ijp]-2*h[ij]+h[ijm])
  }
}

print (c(istep,min(h),mean(h),max(h)))

h2=array(0,c(nx,ny))
for (j in 1:ny){for (i in 1:nx){h2[i,j]=h[i+(j-1)*nx]}}
filled.contour(h2)

}
