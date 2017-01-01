from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import axes3d
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#%reset
#%matplotlib qt

#loadtxt or load...
#dat = loadtxt("data2d.txt")
#bottomtopo = loadtxt("bottomtopo2d.txt")
#par = loadtxt("par2d.txt")

#dat files, like fortran
dat = loadtxt("data2d.dat")
bottomtopo = loadtxt("bottomtopo2d.dat")
par = loadtxt("par2d.dat")


#let's try load
#dat = load("data2d.npy")
#bottomtopo = load("bottomtopo2d.npy")
#par = loadtxt("par2d.txt")

nstop = int(par[0])
dx = par[1]
nx = int(par[2])
Lx = int(par[3])
ny = int(par[4])
Ly = int(par[5])
H0 = int(par[6])


print("--------------PARAMETERS FOR PLOTTING--------------")


X = linspace(0, Lx, nx)
Y = linspace(0, Ly, ny)
X,Y = meshgrid(X,Y)

dat=dat.reshape(nstop,nx,ny)


fig2 = plt.figure(2)
ax2 = fig2.gca(projection='3d')
ax2.set_ylim(0,Ly)
ax2.set_xlim(0,Lx)
ax2.set_xlabel('x [km]')
ax2.set_ylabel('y [km]')
ax2.set_title('Bottom topography')
plt.xticks(linspace(0, Lx, 4),(1/(4*1000))*np.array([0,Lx,Lx*2,Lx*3],dtype=int))
plt.yticks(linspace(0, Ly, 4),(1/(4*1000))*np.array([0,Ly,Ly*2,Ly*3],dtype=int))
ax2.plot_surface(X, Y, bottomtopo, cmap = 'terrain')#rstride=10, cstride=10)
fig2.savefig('figs/bottomtopog.png', bbox_inches='tight')


fig = plt.figure(1) #If i do pcolor, then no need for 3d projection
ax = fig.gca(projection='3d')
#ax = fig.add_subplot(111,projection='3d')

ax.set_xlim(0,Lx)
ax.set_ylim(0,Ly)
ax.set_title('Shallow Water n = 0')
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')
ax.set_zlim(dat.min(),dat.max())
ax.plot_surface(X, Y, dat[0], rstride=10, cstride=10)

#line, = ax.plot_surface(X, Y, dat[0], rstride=2, cstride=2)

def animate(i): #i increment with 1 each step
	ax.clear()
	ax.set_zlim(dat.min(),dat.max())
	ax.set_title('Shallow Water n = {0:3.0f}'.format(i))
	ax.set_ylabel('y [km]')
	ax.set_xlabel('x [km]')
	ax.set_zlim(dat.min(),dat.max())
	plt.xticks(linspace(0, Lx, 4),(1/(4*1000))*np.array([0,Lx,Lx*2,Lx*3],dtype=int))
	plt.yticks(linspace(0, Ly, 4),(1/(4*1000))*np.array([0,Ly,Ly*2,Ly*3],dtype=int))
	#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%5.1f'))
	Z = dat[i]
	#ax.plot_surface(X, Y, bottomtopo, rstride=2, cstride=2, alpha=0.3)
	line = ax.plot_surface(X, Y, dat[i],rstride=10, cstride=10)
	if i == 16:
		fig.savefig('figs/SW2D.png', bbox_inches='tight')
	
	return line,
anim = animation.FuncAnimation(fig,animate, frames = nstop, interval=500)#,blit=False)

plt.show()

