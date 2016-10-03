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
h = zeros((ny,nx),dtype=datatype) #height field
h2 = zeros((ny,nx),dtype=datatype)
box = zeros((30,30),dtype=datatype)
box = box+3999.8
box[1:-2,1:-2] = 4001

u = zeros((ny,nx),dtype=datatype) #horizontal velocity
u = array(u+u0)

v = zeros((ny,nx),dtype=datatype) #horizontal velocity
v = array(v+u0)


ii = zeros((ny,nx),dtype=int)
jj = zeros((ny,nx),dtype=int)
#ix = zeros((ny-2,nx-2),dtype=int)
#jy = zeros((ny-2,nx-2),dtype=int)

#Initialize the grid point index vector, ii,jj
for i in range(nx): 
  ii[i,:] = int(i)#int(i-)
  jj[:,i] = int(i)#int(i-)
#vectorization indices
#ix = ii[:nx+,:nx+]
#jy = jj[:nx+,:nx+]

#Initial h field
h[jj,ii] = H0+1.0*exp(-(ii-nx/5)**2/(nx/3)-(jj-ny/5)**2/(ny/3)) #Gauss wave

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



plt.hold(True)
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
#ax.plot_surface(x[60:70,60:70], y[60:70,60:70], box, color='black')
ax.axis('off')


#we will add a box, it will have 30,30 shape...
#thta means, h(30,30) = 3999

#So, you add boundary conditions at the box sides
testa = zeros((10,10))
for j in range(10):
	for i in range(10):
		testa[j,i] = i
print(testa)
print(testa[0:5,0:5])
#So, h[50:80] will be [50..79]
#so, h[49] = h[48]


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
	# for j in range(1,ny-1):
	# 	for i in range(1,nx-1):
	# 		u2[j,i] = u[j,i]\
	# 			+tau*(u[j,i]*(u[j,i+1]-u[j,i-1])/(2*dx)\
	# 			+v[j,i]*(u[j+1,i]-u[j-1,i])/(2*dy))\
	# 			-tau*g*(h[j,i+1]-h[j,i-1])/(2*dx)

	#127th point, 
	#126 index
	u2[1:ny-1,1:nx-1] = u[1:ny-1,1:nx-1]\
				+tau*(u[1:ny-1,1:nx-1]*(u[1:ny-1,2:nx]-u[1:ny-1,0:nx-2])/(2*dx)\
				+v[1:ny-1,1:nx-1]*(u[2:ny,1:nx-1]-u[0:ny-2,1:nx-1])/(2*dy))\
				-tau*g*(h[1:ny-1,2:nx]-h[1:ny-1,0:nx-2])/(2*dx)


	#You could make v backward too like h, maybe? so only u forward?
	# for j in range(1,ny-1):
	# 	for i in range(1,nx-1):
	# 		v2[j,i] = v[j,i]-tau*(u[j,i]*(v[j,i+1]-v[j,i-1])/(2*dx) \
	# 				+v[j,i]*(v[j+1,i]-v[j-1,i])/(2*dy)) \
	# 				-tau*g*(h[j+1,i]-h[j-1,i])/(2*dy)
	v2[1:ny-1,1:nx-1] = v[1:ny-1,1:nx-1]\
	        -tau*(u[1:ny-1,1:nx-1]*(v[1:ny-1,2:nx]-v[1:ny-1,0:nx-2])/(2*dx) \
			+v[1:ny-1,1:nx-1]*(v[2:ny,1:nx-1]-v[0:ny-2,1:nx-1])/(2*dy)) \
			-tau*g*(h[2:ny,1:nx-1]-h[0:ny-2,1:nx-1])/(2*dy)

	#put boundary conditions for velocity here, they are needed for h2...
	u2[0,:] = 0
	u2[-1,:] = 0
	u2[:,0] = 0
	u2[:,-1] = 0

	v2[0,:] = 0
	v2[-1,:] = 0
	v2[:,0] = 0
	v2[:,-1] = 0

	#Box boundary conditions
	u2[49:82,49:82] = 0
	v2[49:82,49:82] = 0
	#Should include, 49, and 81, right? To go to the actual water edges
	#there, 49 to 82


	#backward, newly updated velocities
	# for j in range(1,ny-1):
	# 	for i in range(1,nx-1):
	# 		h2[j,i] = h[j,i]-tau*u2[j,i]*(h[j,i+1]-h[j,i-1])/(2.0*dx)\
	# 				-tau*v2[j,i]*(h[j+1,i]-h[j-1,i])/(2.0*dy)\
	# 				-tau*h[j,i]*((u2[j,i+1]-u2[j,i-1])/(2.0*dx)+(v2[j+1,i]-v2[j-1,i])/(2.0*dy))

	h2[1:ny-1,1:nx-1] = h[1:ny-1,1:nx-1]\
	                -tau*u2[1:ny-1,1:nx-1]*(h[1:ny-1,2:nx]-h[1:ny-1,0:nx-2])/(2.0*dx)\
					-tau*v2[1:ny-1,1:nx-1]*(h[2:ny,1:nx-1]-h[0:ny-2,1:nx-1])/(2.0*dy)\
					-tau*h[1:ny-1,1:nx-1]*((u2[1:ny-1,2:nx]-u2[1:ny-1,0:nx-2])/(2.0*dx)\
					+(v2[2:ny,1:nx-1]-v2[0:ny-2,1:nx-1])/(2.0*dy))
	
	#boundary conditions
	h2[:,-1] = h2[:,-2] #zero pressure differential at right boundary	
	h2[0,:] = h2[1,:] # dp/dy = 0, at y = 0, using forward differences			
	h2[:,0] = h2[:,1]  #zero pressure differential at left boundary
	h2[-1,:] = h2[-2,:] #zero pressure on that boundary, water does not push into it??

	#box boundary conditions
	h2[50:80,50:80] = 4000 #box

	h2[49:82,49] = h2[49:82,48]
	h2[49:82,81] = h2[49:82,82]	
	h2[49,49:82] = h2[48,49:82]			
	h2[81,49:82] = h2[82,49:82]



	ax.clear()
	ax.axis('off')
	#ax.set_title('Shallow Water n = {0:3.0f}'.format(i))
	#ax.set_ylabel('y [m]')
	#ax.set_xlabel('x [m]')
	ax.set_zlim(3999,4001)
	#ax.plot_surface(x[60:70,60:70], y[60:70,60:70], box, color='black')
	line2 = ax.plot_surface(x[50:80,50:80], y[50:80,50:80], box, rstride=3, cstride=3, color='black')
	line = ax.plot_surface(x, y, h2, rstride=10, cstride=10)
	
	h = h2
	u = u2
	v = v2

	return line2, line

anim = animation.FuncAnimation(fig, animate, frames = nstop, interval=0.5)#,blit=False)
plt.show()

