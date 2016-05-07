import random
import math
import numpy as np
import matplotlib as ml
import matplotlib.pyplot as plt

def find_stack(ij,donor,ndon,stack,nstack):
	for k in range(ndon[ij]):
		ijk=donor[k][ij]
		nstack[0]=nstack[0]+1
		stack[nstack[0]]=ijk
		find_stack (ijk,donor,ndon,stack,nstack)

xl=100.e3
yl=100.e3
nx=101
ny=101
nn=nx*ny
dx=xl/(nx-1)
dy=yl/(ny-1)
dt=1000.
nstep=1000
nfreq=10

h=[]

for i in range(0,nn):
	h.append(random.random())

k=2.e-6
n=2.
m=0.8
u=2.e-3
tol=1.e-3
kd=0.1

for istep in range(0,nstep):

	rec=[i for i in range(nn)]
	length=[0. for i in range(nn)]

	for j in range(1,ny-1):
		for i in range(1,nx-1):
			ij=i+j*nx
			smax=1.e-10
			for jj in range(-1,2):
				for ii in range(-1,2):
					iii=i+ii
					jjj=j+jj
					ijk=iii+jjj*nx
					if ijk != ij:
						l=math.sqrt((dx*ii)**2+(dy*jj)**2)
						slope=(h[ij]-h[ijk])/l
						if slope > smax:
							smax=slope
							rec[ij]=ijk
							length[ij]=l

	ndon=[0 for i in range(nn)]
	donor=[[0 for i in range(nn)] for i in range(8)]
	for ij in range(nn):
		if rec[ij] != ij:
			ijk=rec[ij]
			donor[ndon[ijk]][ijk]=ij
                        ndon[ijk]=ndon[ijk]+1

	stack=[0 for i in range(nn)]
        nstack=[0]
	for ij in range(nn):
		if rec[ij] == ij:
			stack[nstack[0]]=ij
			nstack[0]=nstack[0]+1
			find_stack(ij,donor,ndon,stack,nstack)

	a=[dx*dy for i in range(nn)]
	for ij in range(nn-1,-1,-1):
		ijk=stack[ij]
		if rec[ijk] != ijk:
			a[rec[ijk]]=a[rec[ijk]]+a[ijk]

	for j in range(1,ny-1):
		for i in range(1,nx-1):
			ij=i+j*nx
			h[ij]=h[ij]+u*dt

	for ij in range(nn):
		ijk=stack[ij]
		ijr=rec[ijk]
		if ijr != ijk:
			fact=k*dt*a[ijk]**m/length[ijk]**n
			h0=h[ijk]
			hp=h0
			diff=tol*2
			while abs(diff) > tol:
				h[ijk]=h[ijk]-(h[ijk]-h0+fact*(h[ijk]-h[ijr])**n)/(1.+fact*n*(h[ijk]-h[ijr])**(n-1))
				diff=h[ijk]-hp
				hp=h[ijk]

	factx=kd*dt/dx**2
	facty=kd*dt/dy**2
	for j in range(1,ny-1)
		for i in range(1,nx-1)
			ij=i+j*nx
			imj=i-1+j*nx
			ipj=i+1+j*nx
			ijm=i+(j-1)*nx
			ijp=i+(j+1)*nx
			h[ij]=h[ij]+factx*(h[ipj]-2*h[ij]+h[imj])+facty*(h[ijp]-2*h[ij]+h[ijm])

	if (istep//nfreq)*nfreq == istep:
		print min(h),sum(h)/nn,max(h)
