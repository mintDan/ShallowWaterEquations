from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import axes3d
#import numpy as 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
#%reset
#%matplotlib qt

		
nx = 128											#gridpoints x
ny = nx													#gridpoints y										#gridpoints to make grid periodic
xmin = 0.0
xmax = 1000000.0
Lx = xmax-xmin
dx = Lx/(nx-1)
ymin = xmin
ymax = xmax
Ly = ymax-ymin
dy = Ly/(ny-1)
H0 = 4000.0
nstop = 1000
u0 = 0 #30
v0 = u0
g = 9.82
cg = sqrt(g*H0)
factor = 0.2
tau = factor*(1.0/sqrt(2.0))*dx/(sqrt(u0**2+v0**2)+cg)
Ns = 10


x = linspace(0,Lx,nx)
y = linspace(0,Ly,ny)
x,y = meshgrid(x,y)

datatype=float64



#h is from 0->19
h = zeros((nx,nx),dtype=datatype) #height field
h2 = zeros((nx,nx),dtype=datatype)

u = zeros((nx,nx),dtype=datatype) #horizontal velocity
u = array(u+u0)

v = zeros((nx,nx),dtype=datatype) #horizontal velocity
v = array(v+u0)


ii = zeros((nx,nx),dtype=int)
jj = zeros((nx,nx),dtype=int)


#Initialize the grid point index vector, ii,jj
for i in range(nx): 
  ii[i,:] = int(i)#int(i-)
  jj[:,i] = int(i)#int(i-)
#vectorization indices
#ii = ii[:nx+,:nx+]
#jj = jj[:nx+,:nx+]

#Initial h field
h[ii,jj] = H0+1.0*exp(-(ii-nx/5)**2/(nx/3)-(jj-ny/5)**2/(ny/3)) #Gauss wave

#boundary conditions
u[0,:] = 0
u[:,0] = 0
u[-1,:] = 0
u[:,-1] = 0

v[0,:] = 0
v[:,0] = 0
v[-1,:] = 0
v[:,-1] = 0

h[:,-1] = h[:,-2] 	
h[0,:] = h[1,:]
h[:,0] = h[:,1]
h[-1,:] = h[-2,:]

h2 = h
u2 = u
v2 = v


fig = plt.figure() #If i do pcolor, then no need for 3d projection
ax = fig.gca(projection='3d')
#ax = fig.add_subplot(111,projection='3d')

#ax.set_xlim(0,Lx)
#ax.set_ylim(0,Ly)
#ax.set_title('Shallow Water n = 0')
#ax.set_ylabel('y [m]')
#ax.set_xlabel('x [m]')
ax.set_zlim(3999,4001)
ax.plot_surface(x, y, h, rstride=10, cstride=10)
ax.axis('off')


def animate(i): #i increment with 1 each step

	"""The integration is done in the animate function.
	Uses global variables and arrays, would love to change it"""

	global u, v, h, h2, u2, v2 

	maxvalu = absolute(u).max()
	maxvalv = absolute(v).max()
	maxvalcg = sqrt(g*h.max())
	tau = factor*(1.0/sqrt(2.0))*dx/(sqrt(maxvalu**2+maxvalv**2)+maxvalcg)




	#forward velocity
	#backward difference
	#rane will include nx-2, which is 126...
	#127 is end point, so we can go up to 126
	#we can also calculate at second index, that is 1, since we start at 0
	#so, range (1,nx-2)
	for i in range(1,nx-1):
		for j in range(1,nx-1):
			u2[i,j] = u[i,j]-tau*(u[i,j]*(u[i+1,j]-u[i-1,j])/(2*dx) \
						+v[i,j]*(u[i,j+1]-u[i,j-1])/(2*dy)) \
						-tau*g*(h[i+1,j]-h[i-1,j])/(2*dx)

	#You could make v backward too like h, maybe? so only u forward?
	for i in range(1,nx-1):
		for j in range(1,ny-1):
			v2[i,j] = v[i,j]-tau*(u[i,j]*(v[i+1,j]-v[i-1,j])/(2*dx) \
					+v[i,j]*(v[i,j+1]-v[i,j-1])/(2*dy)) \
					-tau*g*(h[i,j+1]-h[i,j-1])/(2*dy)

	#put boundary conditions for velocity here, they are needed for h2...
	u2[0,:] = 0
	u2[:,0] = 0
	u2[-1,:] = 0
	u2[:,-1] = 0

	v2[0,:] = 0
	v2[:,0] = 0
	v2[-1,:] = 0
	v2[:,-1] = 0


	#backward, newly updated velocities
	for i in range(1,nx-1):
		for j in range(1,ny-1):
			h2[i,j] = h[i,j]-tau*u2[i,j]*(h[i+1,j]-h[i-1,j])/(2.0*dx) \
					-tau*v2[i,j]*(h[i,j+1]-h[i,j-1])/(2.0*dy)\
					-tau*h[i,j]*((u2[i+1,j]-u2[i-1,j])/(2.0*dx)+(v2[i,j+1]-v2[i,j-1])/(2.0*dy))
	
	#boundary conditions

	#

	h2[:,-1] = h2[:,-2] #zero pressure differential at right boundary	
	h2[0,:] = h2[1,:] # dp/dy = 0, at y = 0, using forward differences			
	h2[:,0] = h2[:,1]  #zero pressure differential at left boundary
	h2[-1,:] = h2[-2,:] #zero pressure on that boundary, water does not push into it??

	ax.clear()
	ax.axis('off')
	#ax.set_title('Shallow Water n = {0:3.0f}'.format(i))
	#ax.set_ylabel('y [m]')
	#ax.set_xlabel('x [m]')
	ax.set_zlim(3999,4001)

	line = ax.plot_surface(x, y, h2, rstride=10, cstride=10)

	h = h2
	u = u2
	v = v2
	return line

anim = animation.FuncAnimation(fig, animate, frames = nstop, interval=1)#,blit=False)
plt.show()

