from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import axes3d
#import numpy as 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
#%reset
#%matplotlib qt


#http://stackoverflow.com/questions/1557571/how-to-get-time-of-a-python-program-execution/1557906#1557906
#on windows use time.clock() instead of time.time()
#https://docs.python.org/3/library/time.html
#Men der står det er deprecated, så bruge
#http://stackoverflow.com/questions/7370801/measure-time-elapsed-in-python
#yea, skal vidst bruge process_time
start_time = time.process_time()
#http://stackoverflow.com/questions/673523/how-to-measure-execution-time-of-command-in-windows-command-line
#alternatively, i windows powershell, Measure-Command {ipython SW.py}




def haloval(arr):
		#Sides		
		arr[0:nhalo,:] = arr[nx:nx+nhalo,:]
		arr[nx+nhalo:nx+2*nhalo,:] = arr[nhalo:2*nhalo,:]
	
		arr[:,0:nhalo] = arr[:,nx:nx+nhalo]
		arr[:,nx+nhalo:nx+2*nhalo] = arr[:,nhalo:2*nhalo]
	
		#Corners
		arr[0:nhalo,0:nhalo] = arr[nx:nx+nhalo,nx:nx+nhalo]
		arr[0:nhalo,nx+nhalo:nx+2*nhalo] = arr[nx:nx+nhalo,nhalo:2*nhalo]

		arr[nx+nhalo:nx+2*nhalo,nx+nhalo:nx+2*nhalo] = arr[nhalo:2*nhalo,nhalo:2*nhalo]
		arr[nx+nhalo:nx+2*nhalo,0:nhalo] = arr[nhalo:2*nhalo,nx:nx+nhalo]
		
		return 
		
nx = 128												#gridpoints x
ny = nx													#gridpoints y
nhalo = 3												#gridpoints to make grid periodic
xmin = 0.0
xmax = 1000000.0
Lx = xmax-xmin
dx = Lx/(nx-1)
ymin = xmin
ymax = xmax
Ly = ymax-ymin
dy = Ly/(ny-1)
H0 = 4000.0
nstop = 600
u0 = 2 #30
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
h = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype) #height field
h2 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)

u = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype) #horizontal velocity
u[:,:] = array(u[:,:]+u0)
u2 = u
v = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype) #horizontal velocity
v[:,:] = array(v[:,:]+u0)
v2 = v

ii = zeros((nx+2*nhalo,nx+2*nhalo),dtype=int)
jj = zeros((nx+2*nhalo,nx+2*nhalo),dtype=int)


#Initialize the grid point index vector, ii,jj
for i in range(nhalo,nx+nhalo): 
  ii[i,:] = int(i)#int(i-nhalo)
  jj[:,i] = int(i)#int(i-nhalo)
haloval(ii)
haloval(jj)
#vectorization indices
ix = ii[nhalo:nx+nhalo,nhalo:nx+nhalo]
jy = jj[nhalo:nx+nhalo,nhalo:nx+nhalo]

#Initial h field
h[ix,jy] = H0+1.0*exp(-(ix-nx/2)**2/(nx/3)-(jy-ny/2)**2/(ny/3)) #Gauss wave
haloval(h[:,:])




print('runtime = {0:4.2f}'.format(time.process_time()-start_time))


fig = plt.figure() #If i do pcolor, then no need for 3d projection
ax = fig.gca(projection='3d')
#ax = fig.add_subplot(111,projection='3d')

ax.set_xlim(0,Lx)
ax.set_ylim(0,Ly)
ax.set_title('Shallow Water n = 0')
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')
ax.set_zlim(3999,4001)
ax.plot_surface(x, y, h[nhalo:nx+nhalo,nhalo:nx+nhalo], rstride=10, cstride=10)





def animate(i): #i increment with 1 each step
	global u, v, h, h2, u2, v2 
	#forward velocity
	u2[ix,jy] = u[ix,jy]-tau*(u[ix,jy]*(u[ix+1,jy]-u[ix-1,jy])/(2.0*dx) \
					+v[ix,jy]*(u[ix,jy+1]-u[ix,jy-1])/(2.0*dy)) \
					-tau*g*(h[ix+1,jy]-h[ix-1,jy])/(2.0*dx)
	#You could make v backward too like h, maybe? so only u forward?
	v2[ix,jy] = v[ix,jy]-tau*(u[ix,jy]*(v[ix+1,jy]-v[ix-1,jy])/(2.0*dx) \
					+v[ix,jy]*(v[ix,jy+1]-v[ix,jy-1])/(2.0*dy)) \
					-tau*g*(h[ix,jy+1]-h[ix,jy-1])/(2.0*dy)
	
	#backward, newly updated velocities
	h2[ix,jy] = h[ix,jy]-tau*u[ix,jy]*(h[ix+1,jy]-h[ix-1,jy])/(2.0*dx) \
					-tau*v[ix,jy]*(h[ix,jy+1]-h[ix,jy-1])/(2.0*dy)\
					-tau*h[ix,jy]*((u[ix+1,jy]-u[ix-1,jy])/(2.0*dx)+(v[ix,jy+1]-v[ix,jy-1])/(2.0*dy))
	
	haloval(u2[:,:])
	haloval(v2[:,:])
	haloval(h2[:,:])

	ax.clear()
	ax.set_title('Shallow Water n = {0:3.0f}'.format(i))
	ax.set_ylabel('y [m]')
	ax.set_xlabel('x [m]')
	ax.set_zlim(3999,4001)

	line = ax.plot_surface(x, y, h2[nhalo:nx+nhalo,nhalo:nx+nhalo], rstride=10, cstride=10)

	h = h2
	u = u2
	v = v2
	return line

anim = animation.FuncAnimation(fig, animate, frames = nstop, interval=4)#,blit=False)
plt.show()

