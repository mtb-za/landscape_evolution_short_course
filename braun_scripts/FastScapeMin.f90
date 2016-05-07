! test program for FastScape Course

program FastScape0

implicit none

! declaring arrays

real, dimension(:), allocatable :: h,a,length
real, dimension(:,:), allocatable :: x,y,z,profile
integer, dimension(:), allocatable :: rec,ndon,stack
integer, dimension(:,:), allocatable :: donor

integer nx,ny,nn,nstep,nfreq,nstack
integer i,j,ij,ii,jj,iii,jjj,ijk,ijr,istep,iout,ip,ipmax
integer imj,ipj,ijm,ijp
real xl,yl,dx,dy,dt,k,n,m,u,l,slope,smax,kd,precipitation
real diff,fact,h0,hp,tol,factx,facty

! defining size of the problem

nx=51;ny=51
nn=nx*ny

! allocating memory

allocate (h(nn),a(nn),length(nn),rec(nn),ndon(nn),stack(nn),donor(8,nn))
allocate (x(nx,ny),y(nx,ny),z(nx,ny))

! defining geometrical and temporal constants

xl=100.e3;yl=100.e3
dx=xl/(nx-1);dy=yl/(ny-1)
dt=10000.
nstep=100
nfreq=10
allocate (profile(nn,nstep/nfreq))
profile=0.
iout=0
ipmax=0

! generating initial topography

call random_number (h)
  do j=1,ny
    do i=1,nx
    ij=i+(j-1)*nx
    if (i.eq.1.or.i.eq.nx.or.j.eq.1.or.j.eq.ny) h(ij)=0.d0
    enddo
  enddo

! initializing erosional parameters

k=1.e-4;n=1;m=0.4
u=2.e-3
tol=1.e-3
kd=0.
precipitation=1.

! begining of time stepping

  do istep=1,nstep

! initializing rec and length

    do ij=1,nn
    rec(ij)=ij
    length(ij)=0.
    enddo

! computing receiver array

    do j=2,ny-1
      do i=2,nx-1
      ij=i+(j-1)*nx
      smax=tiny(smax)
        do jj=-1,1
          do ii=-1,1
          iii=i+ii
          jjj=j+jj
          ijk=iii+(jjj-1)*nx
            if (ijk.ne.ij) then
            l=sqrt((dx*ii)**2+(dy*jj)**2)
            slope=(h(ij)-h(ijk))/l
              if (slope.gt.smax) then
              smax=slope
              rec(ij)=ijk
              length(ij)=l
              endif
            endif
          enddo
        enddo
      enddo
    enddo

! initialising number of donors per node to 0

  ndon=0

! computing donor arrays

    do ij=1,nn
      if (rec(ij).ne.ij) then
      ijk=rec(ij)
      ndon(ijk)=ndon(ijk)+1
      donor(ndon(ijk),ijk)=ij
      endif
    enddo

! computing stack

  nstack=0
    do ij=1,nn
      if (rec(ij).eq.ij) then
      nstack=nstack+1
      stack(nstack)=ij
      call find_stack (ij,donor,ndon,nn,stack,nstack)
      endif
    enddo

! computing drainage area

  a=dx*dy*precipitation
    do ij=nn,1,-1
    ijk=stack(ij)
      if (rec(ijk).ne.ijk) then
      a(rec(ijk))=a(rec(ijk))+a(ijk)
      endif
    enddo

! adding uplift to landscape

    do j=2,ny-1
      do i=2,nx-1
      ij=i+(j-1)*nx
      h(ij)=h(ij)+u*dt
      enddo
    enddo

! computing erosion

    do ij=1,nn
    ijk=stack(ij)
    ijr=rec(ijk)
      if (ijr.ne.ijk) then
      fact=k*dt*a(ijk)**m/length(ijk)**n
      h(ijk)=(h(ijk)+fact*h(ijr))/(1.+fact)
      endif
    enddo

! adding diffusion

factx=kd*dt/dx**2
facty=kd*dt/dy**2
  do j=2,ny-1
    do i=2,nx-1
    ij=i+(j-1)*nx
    imj=i-1+(j-1)*nx
    ipj=i+1+(j-1)*nx
    ijm=i+(j-2)*nx
    ijp=i+j*nx
    h(ij)=h(ij)+factx*(h(ipj)-2.*h(ij)+h(imj))+facty*(h(ijp)-2.*h(ij)+h(ijm))
    enddo
  enddo

! outputing summary of results every nfreq steps

    if ((istep/nfreq)*nfreq.eq.istep) then
    print*,minval(h),sum(h)/nn,maxval(h)
    iout=iout+1
    i=maxloc(a,1)
    ip=0
      do while (ndon(i).ne.0)
      i=donor(maxloc(a(donor(1:ndon(i),i)),1),i)
      ip=ip+1
      profile(ip,iout)=h(i)
      enddo
    ipmax=max(ipmax,ip)
    endif

  enddo

! output final topography to ascii file

open (7,file='FinalTopo.txt',status='unknown')
  do j=1,ny
  write (7,*) (h(i+(j-1)*nx),i=1,nx)
  enddo
close (7)

! output largest river profile to ascii file

print*,ipmax,iout
open (7,file='Profile.txt',status='unknown')
  do i=1,ipmax
  write (7,*) (profile(i,j),j=1,iout)
  enddo
close (7)

end

!----------

! recursive routine to compute the stack

recursive subroutine find_stack	(ij,donor,ndon,nn,stack,nstack)

implicit none

integer ij,k,ijk,nn,nstack
integer donor(8,nn),ndon(nn),stack(nstack)

  do k=1,ndon(ij)
  ijk=donor(k,ij)
  nstack=nstack+1
  stack(nstack)=ijk
  call find_stack (ijk,donor,ndon,nn,stack,nstack)
  enddo

return
end
