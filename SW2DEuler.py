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



def step(h,u,v):
	#updating short timestep tau, based on u,v and h.
	#Updating tau, of course if u,v becomes negative, then we want .min()
	#so need to fix this, for peculiar cases.
	maxvalu = absolute(u[nhalo:nx+nhalo,nhalo:nx+nhalo]).max()
	maxvalv = absolute(v[nhalo:nx+nhalo,nhalo:nx+nhalo]).max()
	maxvalcg = sqrt(g*h[nhalo:nx+nhalo,nhalo:nx+nhalo].max())
	tau = factor*(1.0/sqrt(2.0))*dx/(sqrt(maxvalu**2+maxvalv**2)+cg)
	
	
	#Should see if matrix multiplication is faster than all this?
	#And maybe write u,v better, and maybe, make a single .array...
	#forward velocity
	u[ix,jy] = u[ix,jy]-tau*(u[ix,jy]*(u[ix+1,jy]-u[ix-1,jy])/(2.0*dx) \
					+v[ix,jy]*(u[ix,jy+1]-u[ix,jy-1])/(2.0*dy)) \
					-tau*g*(h[ix+1,jy]-h[ix-1,jy])/(2.0*dx)
	#You could make v backward too like h, maybe? so only u forward?
	v[ix,jy] = v[ix,jy]-tau*(u[ix,jy]*(v[ix+1,jy]-v[ix-1,jy])/(2.0*dx) \
					+v[ix,jy]*(v[ix,jy+1]-v[ix,jy-1])/(2.0*dy)) \
					-tau*g*(h[ix,jy+1]-h[ix,jy-1])/(2.0*dy)
	
	#backward, newly updated velocities
	h[ix,jy] = h[ix,jy]-tau*u[ix,jy]*(h[ix+1,jy]-h[ix-1,jy])/(2.0*dx) \
					-tau*v[ix,jy]*(h[ix,jy+1]-h[ix,jy-1])/(2.0*dy)\
					-tau*h[ix,jy]*((u[ix+1,jy]-u[ix-1,jy])/(2.0*dx)+(v[ix,jy+1]-v[ix,jy-1])/(2.0*dy))
	#Should change h to d, but maybe not all places...
	
	haloval(u[:,:])
	haloval(v[:,:])
	haloval(h[:,:])

	return h,u,v

print('runtime = {0:4.2f}'.format(time.process_time()-start_time))


fig = plt.figure() #If i do pcolor, then no need for 3d projection
ax = fig.gca(projection='3d')
#ax = fig.add_subplot(111,projection='3d')

ax.set_xlim(0,Lx)
ax.set_ylim(0,Ly)
ax.set_title('Shallow Water n = 0')
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')
ax.set_zlim(h.min(),h.max())
#ax.xaxis.set_ticks(0,Lx,4)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%5.1f'))
ax.plot_surface(x, y, h[nhalo:nx+nhalo,nhalo:nx+nhalo], rstride=10, cstride=10)

#line, = ax.plot_surface(X, Y, dat[0], rstride=2, cstride=2)



#contour ,cmap='jet'.... 
#pcolor looks like matlab equilevant contour
#http://matplotlib.org/examples/images_contours_and_fields/contourf_log.html
#This also looks close to Matlab
#http://stackoverflow.com/questions/15601096/contour-graph-in-python
#http://matplotlib.org/examples/pylab_examples/contourf_demo.html

# plt.pcolor(X,Y,dat[0],cmap='jet')
# plt.hold(False)
# plt.xlabel('x [m]')
# plt.ylabel('y [m]')
# plt.colorbar()
def animate(h,u,v,h2,u2,v2): #i increment with 1 each step
	h2,u2,v2 = step(h,u,v)
	ax.clear()
	#ax.set_zlim(h[i,nhalo:nx+nhalo,nhalo:nx+nhalo].min(),dat.max())
	ax.set_title('Shallow Water n = {0:3.0f}'.format(i))
	ax.set_ylabel('y [m]')
	ax.set_xlabel('x [m]')
	ax.set_zlim(h.min(),h.max())
	#ax.plot_surface(X, Y, bottomtopo, rstride=2, cstride=2, alpha=0.3)
	line = ax.plot_surface(x, y, h[nhalo:nx+nhalo,nhalo:nx+nhalo], rstride=10, cstride=10)
	#plot=plt.contour(X,Y,dat[i],3,cmap='jet')
	#plt.pcolor(X,Y,dat[i],cmap='jet')
	h = h2
	u = u2
	v = v2
	return line,

anim = animation.FuncAnimation(fig,animate(h,u,v,h2,u2,v2), frames = nstop, interval=4)#,blit=False)
plt.show()

