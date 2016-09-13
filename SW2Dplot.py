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


fig = plt.figure() #If i do pcolor, then no need for 3d projection
ax = fig.gca(projection='3d')
#ax = fig.add_subplot(111,projection='3d')

ax.set_xlim(0,Lx)
ax.set_ylim(0,Ly)
ax.set_title('Shallow Water n = 0')
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')
ax.set_zlim(dat.min(),dat.max())
#ax.xaxis.set_ticks(0,Lx,4)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%5.1f'))
ax.plot_surface(X, Y, dat[0], rstride=10, cstride=10)

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
def animate(i): #i increment with 1 each step
	ax.clear()
	ax.set_zlim(dat.min(),dat.max())
	ax.set_title('Shallow Water n = {0:3.0f}'.format(i))
	ax.set_ylabel('y [m]')
	ax.set_xlabel('x [m]')
	ax.set_zlim(dat.min(),dat.max())
	#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%5.1f'))
	Z = dat[i]
	#ax.plot_surface(X, Y, bottomtopo, rstride=2, cstride=2, alpha=0.3)
	line = ax.plot_surface(X, Y, dat[i], rstride=10, cstride=10)
	#plot=plt.contour(X,Y,dat[i],3,cmap='jet')
	#plt.pcolor(X,Y,dat[i],cmap='jet')
	return line,
anim = animation.FuncAnimation(fig,animate, frames = nstop, interval=500)#,blit=False)

plt.show()

#figure(1)
#for n in range(nt-1):
	#plt.plot(x,h[n,nhalo:nx+nhalo],'blue')
	#plt.show()
#figure(1)
#plt.ion()
# for n in range(nt-1):
	#ax.clear()
	# ax.clear()
	# plt.plot(x,y+n)
	# plt.draw()
	# plt.show()
	#plt.close()
#plt.show()
