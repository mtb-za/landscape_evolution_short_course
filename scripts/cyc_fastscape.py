from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

debug = False
stop_uplift = False
nx,ny = 101, 101
nn = nx*ny
xl, yl = 100000,100000
kf = 0.0001
kd = 10 #diffusion constant
n = 1
m = 0.4
u = 0.002
dt = 10000
dx = xl/(nx-1) #change in x
dy = yl/(ny-1) #change in y

def find_stack(ij, stack, nstack, donors, num_donors, debug):
    '''
    Recursively put the donors of ij on the stack.
    '''
    for k in range(num_donors[ij]):
        if debug:
            print("ij:", ij, "k:", k, "nstack:", nstack)
        ijk = donors[k][ij]
        stack[nstack[0]] = ijk
        nstack[0] = nstack[0] + 1
        if debug:
            print ("Before: ",ij, "\t", k, "\t", nstack, "\n", stack) #Debugging comment
        find_stack(ijk, stack, nstack, donors, num_donors, debug)

np.random.seed(123) #For checking that things are running the same each time.
h = []
for i in range(nn):
    h.append(np.random.random())

for outer_loop in range(101):
	if outer_loop == 50 and stop_uplift:
		u = 0
	
	rec_array = [i for i in range(nn)]
	length =  [0 for i in range(nn)]
	for i in range(nx):
		for j in range(1, ny-1):
			ij = i + j * nx
			max_slope = 0
			for ii in range(-1, 2):
				for jj in range(-1, 2):
					ii_val = i + ii
					jj_val = j + jj
					if ii_val < 0:
						ii_val = nx-1
					elif ii_val > nx-1:
						ii_val = 0
					ijk = ii_val + jj_val * nx
					if ijk != ij:
						distance = np.sqrt(((dx**2)*(ii**2))+((dy**2)*(jj**2)))
						slope = (h[ij] - h[ijk]) / distance
						if slope > max_slope:
							max_slope = slope
							rec_array[ij] = ijk
							length[ij] = distance

	num_donors = [0 for i in range(nn)]
	donors = [[0 for i in range(nn)] for j in range(8)]

	for ij in range(nn):
		if rec_array[ij] != ij:
			k = rec_array[ij]
			donors[num_donors[k]][k] = ij
			num_donors[k] += 1

	#stack
	nstack = [0]
	stack = [0 for i in range(nn)]
	for ij in range(nn):
		if rec_array[ij] == ij:
			stack[nstack[0]] = ij
			nstack[0] = nstack[0] + 1
			if debug:
				print("ij:", ij, "nstack:", nstack) #debugging comment
			find_stack(ij, stack, nstack, donors, num_donors, debug)

	rev_stack = []
	for i in range(nn):
		if debug:
			print(nn-i-1, stack[nn-i-1])
		rev_stack.append(stack[nn-i-1])

	#area
	area = [dx*dy for i in range(nn)]
	for ij in rev_stack:
		if ij != rec_array[ij]:
			area[rec_array[ij]] += area[ij]

	#uplift
	for i in range(nx):
		for j in range(1, ny-1):
			ij = i+j*nx
			h[ij] = h[ij] + u * dt
		
	#erode
	for ij in range(nn):
		ijk = stack[ij]
		if rec_array[ijk] != ijk:
			f = area[ij]**m*kf*dt / length[ijk]
			h[ijk] = ( h[ijk] + f*h[rec_array[ijk]] ) / (1+f)
	
	
	#stuff for diffusion:
	hp = []
	hp = [h[ij] for ij in range(nn)]
	
	fact_x = dt*kd/dx**2
	fact_y = dt*kd/dy**2
	
	for j in range(1,ny-1):
		for i in range(1,nx-1):
			ipj = i+j*nx
			imj = -i+j*nx
			ijp = i+(j+1)*nx
			if ijp > ny-1:
				ijp = 0
			elif ijp < 0:
				ijp = ny-1
			ijm = i+(j-1)*nx
			if ijm < 0:
				ijm = ny-1
			
			h[ij] = h[ij] + fact_x * ( hp[ipj] - 2*hp[ij] + hp[ijm] ) + fact_y * ( hp[ijp] - 2*hp[ij] + hp[ijm] )
	
	if outer_loop % 10 == 0:
		print("snapshotting at", outer_loop, "iterations")
		to_plot = np.array(h)
		to_plot = np.reshape(h, (nx,ny))
		x = range(nx)
		y = range(ny)

		X, Y = np.meshgrid(x, y)

		plt.contourf(X, Y, to_plot, cmap='viridis')
		filename = str(outer_loop)+".png"
		plt.savefig(filename)
		if outer_loop == 100:
			fig = plt.figure()
			ax = fig.add_subplot(1,1,1, projection='3d')
			surf = ax.plot_surface(X, Y, to_plot, rstride=1, cstride=1, linewidth=0, cmap='viridis')
			plt.show()
