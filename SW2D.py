from numpy import *
from pylab import *
#from decimal import Decimal
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
#%reset
#%matplotlib qt

#Checking for runtime
start_time = time.process_time()

#Normally working, from fortran
f1 = open('par2d.dat','a')
f2 = open('data2d.dat','a')
f3 = open('bottomtopo2d.dat','a')
f4 = open('totalmass.dat','a')


#Trying to do it without for loops, numpy.save or numpy.savetxt
#savetxt and loadtxt
#http://stackoverflow.com/questions/27786868/python3-numpy-appending-to-a-file-using-numpy-savetxt
#Should open it in binary mode?
#f = open('par2d.txt','ab')
#f2 = open('data2d.dat','ab') #WORKS
#f3 = open('bottomtopo2d.txt','ab')
#http://stackoverflow.com/questions/24691755/how-to-format-in-numpy-savetxt-such-that-zeros-are-saved-only-as-0
#You can add format to savetxt... but it seems like savetxt source code is itself for loops
#so maybe that's why my own is just as fast..




#save
#http://stackoverflow.com/questions/25749215/save-a-numpy-matrix
#f = open('par2d.npy','a')
#f2 = open('data2d.npy','ab')
#f3 = open('bottomtopo2d.npy','ab')


#https://mail.scipy.org/pipermail/numpy-discussion/2007-November/030099.html
#savetxt again
#f2 = open('data2d.dat','a')


#https://www.reddit.com/r/Python/comments/2szc8r/saving_multiple_2d_matrices_to_a_single_file/
#savetxt igenigen
#f2 = open("my_file.txt", "a")

#call full array h[n2,:] since it's the endpoints we're changing
#TEST this haloval if it gives correct arrays!!!
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



	
def bicubicvec(a,b,p,q,arr,arr2):
	#Maybe, if you ONLY do [ix,jy] on a multidimensional array,
	#it looks at all indices i,j as ix,jy... a short cut perhaps?
	#If that's the case, then maybe actually, python arrays are better than fortran arrays.
	
	wa1[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]-1]/6. \
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]-1]/2. \
							+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]-1]/2. \
							-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]-1]/6.
						
	wa2[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]]/6.\
							+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]]/2.\
							+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]]/2.\
							-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]]/6.

	wa3[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+1]/6.\
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+1]/2.\
						+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+1]/2.\
						-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+1]/6.

	wa4[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+2]/6.\
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+2]/2.\
						+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+2]/2.\
						-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+2]/6.

	arr2[ix,jy] = -(b[ix,jy])*(1.-b[ix,jy])*(2.-b[ix,jy])*wa1[ix,jy]/6.\
							+(1.-b[ix,jy]**2)*(2.-b[ix,jy])*wa2[ix,jy]/2.\
							+b[ix,jy]*(1.+b[ix,jy])*(2.-b[ix,jy])*wa3[ix,jy]/2.\
							-b[ix,jy]*(1.-b[ix,jy]**2)*wa4[ix,jy]/6.			

	haloval(arr2)
	return


def bicubicLMCSLvec(a,b,p,q,Cupx,Cupy,arr,arr2):
	#nx, ny ix suppose, or nhalo,nx+nhalo...
	#Is it faster to do in place operation, like mass *=0? Rather than creating the array again?
	mass = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype) 
	#Quasi-Lagrangian areas at the departure points
	a1[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy],q[ix,jy]-1]+Cupx[p[ix,jy]-2,q[ix,jy]-1])**2+(Cupy[p[ix,jy],q[ix,jy]-1]-Cupy[p[ix,jy]-2,q[ix,jy]-1])**2)\
			 *sqrt((2.0*dy - Cupy[p[ix,jy]-1,q[ix,jy]]+Cupy[p[ix,jy]-1,q[ix,jy]-2])**2+(Cupx[p[ix,jy]-1,q[ix,jy]]-Cupx[p[ix,jy]-1,q[ix,jy]-2])**2)
		
	a2[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+1,q[ix,jy]-1]+Cupx[p[ix,jy]-1,q[ix,jy]-1])**2+(Cupy[p[ix,jy]+1,q[ix,jy]-1]-Cupy[p[ix,jy]-1,q[ix,jy]-1])**2)\
			 *sqrt((2.0*dy - Cupy[p[ix,jy],q[ix,jy]]+Cupy[p[ix,jy],q[ix,jy]-2])**2+(Cupx[p[ix,jy]-1,q[ix,jy]]-Cupx[p[ix,jy]-1,q[ix,jy]-2])**2)		

	a3[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+2,q[ix,jy]-1]+Cupx[p[ix,jy],q[ix,jy]-1])**2+(Cupy[p[ix,jy]+2,q[ix,jy]-1]-Cupy[p[ix,jy],q[ix,jy]-1])**2)\
			 *sqrt((2.0*dy - Cupy[p[ix,jy]+1,q[ix,jy]]+Cupy[p[ix,jy]+1,q[ix,jy]-2])**2+(Cupx[p[ix,jy]+1,q[ix,jy]]-Cupx[p[ix,jy]+1,q[ix,jy]-2])**2)	
	
	a4[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+3,q[ix,jy]-1]+Cupx[p[ix,jy]+1,q[ix,jy]-1])**2+(Cupy[p[ix,jy]+3,q[ix,jy]-1]-Cupy[p[ix,jy]+1,q[ix,jy]-1])**2)\
			 *sqrt((2.0*dy - Cupy[p[ix,jy]+2,q[ix,jy]]+Cupy[p[ix,jy]+2,q[ix,jy]-2])**2+(Cupx[p[ix,jy]+2,q[ix,jy]]-Cupx[p[ix,jy]+2,q[ix,jy]-2])**2)	

	a5[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+1,q[ix,jy]]+Cupx[p[ix,jy]-1,q[ix,jy]])**2+(Cupy[p[ix,jy]+1,q[ix,jy]]-Cupy[p[ix,jy]-1,q[ix,jy]])**2)\
			 *sqrt((2.0*dy - Cupy[p[ix,jy]-1,q[ix,jy]+1]+Cupy[p[ix,jy]-1,q[ix,jy]-1])**2+(Cupx[p[ix,jy]-1,q[ix,jy]+1]-Cupx[p[ix,jy]-1,q[ix,jy]-1])**2)	
	
	a6[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+1,q[ix,jy]]+Cupx[p[ix,jy]-1,q[ix,jy]])**2+(Cupy[p[ix,jy]+1,q[ix,jy]]-Cupy[p[ix,jy]-1,q[ix,jy]])**2)\
			 *sqrt((2.0*dy - Cupy[p[ix,jy],q[ix,jy]+1]+Cupy[p[ix,jy],q[ix,jy]-1])**2+(Cupx[p[ix,jy],q[ix,jy]+1]-Cupx[p[ix,jy],q[ix,jy]-1])**2)	

	a7[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+2,q[ix,jy]]+Cupx[p[ix,jy],q[ix,jy]])**2+(Cupy[p[ix,jy]+2,q[ix,jy]]-Cupy[p[ix,jy],q[ix,jy]])**2)\
			 *sqrt((2.0*dy - Cupy[p[ix,jy]+1,q[ix,jy]+1]+Cupy[p[ix,jy]+1,q[ix,jy]-1])**2+(Cupx[p[ix,jy]+1,q[ix,jy]+1]-Cupx[p[ix,jy]+1,q[ix,jy]-1])**2)	
	
	a8[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+3,q[ix,jy]]+Cupx[p[ix,jy]+1,q[ix,jy]])**2+(Cupy[p[ix,jy]+3,q[ix,jy]]-Cupy[p[ix,jy]+1,q[ix,jy]])**2)\
			 *sqrt((2.0*dy - Cupy[p[ix,jy]+2,q[ix,jy]+1]+Cupy[p[ix,jy]+2,q[ix,jy]-1])**2+(Cupx[p[ix,jy]+2,q[ix,jy]+1]-Cupx[p[ix,jy]+2,q[ix,jy]-1])**2)
	
	a9[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy],q[ix,jy]+1]+Cupx[p[ix,jy]-2,q[ix,jy]+1])**2+(Cupy[p[ix,jy],q[ix,jy]+1]-Cupy[p[ix,jy]-2,q[ix,jy]+1])**2)\
			 *sqrt((2.0*dy - Cupy[p[ix,jy]-1,q[ix,jy]+2]+Cupy[p[ix,jy]-1,q[ix,jy]])**2+(Cupx[p[ix,jy]-1,q[ix,jy]+2]-Cupx[p[ix,jy]+2,q[ix,jy]])**2)
	
	a10[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+1,q[ix,jy]+1]+Cupx[p[ix,jy]-1,q[ix,jy]+1])**2+(Cupy[p[ix,jy]+1,q[ix,jy]+1]-Cupy[p[ix,jy]-1,q[ix,jy]+1])**2)\
			  *sqrt((2.0*dy - Cupy[p[ix,jy],q[ix,jy]+2]+Cupy[p[ix,jy],q[ix,jy]])**2+(Cupx[p[ix,jy],q[ix,jy]+2]-Cupx[p[ix,jy],q[ix,jy]])**2)		
	
	a11[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+2,q[ix,jy]+1]+Cupx[p[ix,jy],q[ix,jy]+1])**2+(Cupy[p[ix,jy]+2,q[ix,jy]+1]-Cupy[p[ix,jy],q[ix,jy]+1])**2)\
			  *sqrt((2.0*dy - Cupy[p[ix,jy]+1,q[ix,jy]+2]+Cupy[p[ix,jy]+1,q[ix,jy]])**2+(Cupx[p[ix,jy]+1,q[ix,jy]+2]-Cupx[p[ix,jy]+1,q[ix,jy]])**2)				
	
	a12[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+3,q[ix,jy]+1]+Cupx[p[ix,jy]+1,q[ix,jy]+1])**2+(Cupy[p[ix,jy]+3,q[ix,jy]+1]-Cupy[p[ix,jy]+1,q[ix,jy]+1])**2)\
			  *sqrt((2.0*dy - Cupy[p[ix,jy]+2,q[ix,jy]+2]+Cupy[p[ix,jy]+2,q[ix,jy]])**2+(Cupx[p[ix,jy]+2,q[ix,jy]+2]-Cupx[p[ix,jy]+2,q[ix,jy]])**2)		

	a13[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy],q[ix,jy]+2]+Cupx[p[ix,jy]-2,q[ix,jy]+2])**2+(Cupy[p[ix,jy],q[ix,jy]+2]-Cupy[p[ix,jy]-2,q[ix,jy]+2])**2)\
			  *sqrt((2.0*dy - Cupy[p[ix,jy]-1,q[ix,jy]+3]+Cupy[p[ix,jy]-1,q[ix,jy]+1])**2+(Cupx[p[ix,jy]-1,q[ix,jy]+3]-Cupx[p[ix,jy]-1,q[ix,jy]+1])**2)		

	a14[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+1,q[ix,jy]+2]+Cupx[p[ix,jy]-1,q[ix,jy]+2])**2+(Cupy[p[ix,jy]+1,q[ix,jy]+2]-Cupy[p[ix,jy]-1,q[ix,jy]+2])**2)\
			  *sqrt((2.0*dy - Cupy[p[ix,jy],q[ix,jy]+3]+Cupy[p[ix,jy],q[ix,jy]+1])**2+(Cupx[p[ix,jy],q[ix,jy]+3]-Cupx[p[ix,jy],q[ix,jy]+1])**2)	

	a15[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+2,q[ix,jy]+2]+Cupx[p[ix,jy],q[ix,jy]+2])**2+(Cupy[p[ix,jy]+2,q[ix,jy]+2]-Cupy[p[ix,jy],q[ix,jy]+2])**2)\
			  *sqrt((2.0*dy - Cupy[p[ix,jy]+1,q[ix,jy]+3]+Cupy[p[ix,jy]+1,q[ix,jy]+1])**2+(Cupx[p[ix,jy]+1,q[ix,jy]+3]-Cupx[p[ix,jy]+1,q[ix,jy]+1])**2)			
	
	a16[ix,jy] = 0.25*sqrt((2.0*dx - Cupx[p[ix,jy]+3,q[ix,jy]+2]+Cupx[p[ix,jy]+1,q[ix,jy]+2])**2+(Cupy[p[ix,jy]+3,q[ix,jy]+2]-Cupy[p[ix,jy]+1,q[ix,jy]+2])**2)\
			  *sqrt((2.0*dy - Cupy[p[ix,jy]+2,q[ix,jy]+3]+Cupy[p[ix,jy]+2,q[ix,jy]+1])**2+(Cupx[p[ix,jy]+2,q[ix,jy]+3]-Cupx[p[ix,jy]+2,q[ix,jy]+1])**2)

			  
	#w1,fixrst y
	mass[p[ix,jy]-1,q[ix,jy]-1] = mass[p[ix,jy]-1,q[ix,jy]-1]-((a[ix,jy])*(1.-a[ix,jy])*(2-a[ix,jy])*(a1[ix,jy]/a0)/6.) \
													*(-(b[ix,jy])*(1.-b[ix,jy])*(2.-b[ix,jy])/6.)
	
	mass[p[ix,jy],q[ix,jy]-1] = mass[p[ix,jy],q[ix,jy]-1]+((1-a[ix,jy]**2)*(2-a[ix,jy])*(a2[ix,jy]/a0)/2.)  \
													*(-(b[ix,jy])*(1.-b[ix,jy])*(2.-b[ix,jy])/6.)
								
	mass[p[ix,jy]+1,q[ix,jy]-1] = mass[p[ix,jy]+1,q[ix,jy]-1]+(a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*(a3[ix,jy]/a0)/2.) \
													*(-(b[ix,jy])*(1.-b[ix,jy])*(2.-b[ix,jy])/6.)

	mass[p[ix,jy]+2,q[ix,jy]-1] = mass[p[ix,jy]+2,q[ix,jy]-1]-(a[ix,jy]*(1.-a[ix,jy]**2)*(a4[ix,jy]/a0)/6.) \
													*(-(b[ix,jy])*(1.-b[ix,jy])*(2.-b[ix,jy])/6.)	
										
	#w2, second y										
	mass[p[ix,jy]-1,q[ix,jy]] = mass[p[ix,jy]-1,q[ix,jy]]-(a[ix,jy]*(1.-a[ix,jy])*(2.-a[ix,jy])*(a5[ix,jy]/a0)/6.)\
													*(+(1.-b[ix,jy]**2)*(2.-b[ix,jy])/2.)	

	mass[p[ix,jy],q[ix,jy]] = mass[p[ix,jy],q[ix,jy]]+((1.-a[ix,jy]**2)*(2.-a[ix,jy])*(a6[ix,jy]/a0)/2.)\
													*(+(1.-b[ix,jy]**2)*(2.-b[ix,jy])/2.)			
	
	mass[p[ix,jy]+1,q[ix,jy]] = mass[p[ix,jy]+1,q[ix,jy]]+(a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*(a7[ix,jy]/a0)/2.)\
													*(+(1.-b[ix,jy]**2)*(2.-b[ix,jy])/2.)	
												
	mass[p[ix,jy]+2,q[ix,jy]] = mass[p[ix,jy]+2,q[ix,jy]]-(a[ix,jy]*(1.-a[ix,jy]**2)*(a8[ix,jy]/a0)/6.)\
													*(+(1.-b[ix,jy]**2)*(2.-b[ix,jy])/2.)				

	#w3, thixrd y								
	mass[p[ix,jy]-1,q[ix,jy]+1] = mass[p[ix,jy]-1,q[ix,jy]+1]-(a[ix,jy]*(1.-a[ix,jy])*(2.-a[ix,jy])*(a9[ix,jy]/a0)/6.)\
													*(+b[ix,jy]*(1.+b[ix,jy])*(2.-b[ix,jy])/2.)
												
	mass[p[ix,jy],q[ix,jy]+1] = mass[p[ix,jy]  ,q[ix,jy]+1]+((1.-a[ix,jy]**2)*(2.-a[ix,jy])*(a10[ix,jy]/a0)/2.)\
													*(+b[ix,jy]*(1.+b[ix,jy])*(2.-b[ix,jy])/2.)

	mass[p[ix,jy]+1,q[ix,jy]+1] = mass[p[ix,jy]+1,q[ix,jy]+1]+(a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*(a11[ix,jy]/a0)/2.)\
													*(+b[ix,jy]*(1.+b[ix,jy])*(2.-b[ix,jy])/2.)

	mass[p[ix,jy]+2,q[ix,jy]+1] = mass[p[ix,jy]+2,q[ix,jy]+1]-(a[ix,jy]*(1.-a[ix,jy]**2)*(a12[ix,jy]/a0)/6.)\
													*(+b[ix,jy]*(1.+b[ix,jy])*(2.-b[ix,jy])/2.)

	#w4 fourth
	mass[p[ix,jy]-1,q[ix,jy]+2] = mass[p[ix,jy]-1,q[ix,jy]+2]-(a[ix,jy]*(1.-a[ix,jy])*(2.-a[ix,jy])*(a13[ix,jy]/a0)/6.)\
													*(-b[ix,jy]*(1.-b[ix,jy]**2)/6.)	
												
	mass[p[ix,jy],q[ix,jy]+2] = mass[p[ix,jy],q[ix,jy]+2]+((1.-a[ix,jy]**2)*(2.-a[ix,jy])*(a14[ix,jy]/a0)/2.)\
													*(-b[ix,jy]*(1.-b[ix,jy]**2)/6.)
	
	mass[p[ix,jy]+1,q[ix,jy]+2] = mass[p[ix,jy]+1,q[ix,jy]+2]+(a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*(a15[ix,jy]/a0)/2.)\
													*(-b[ix,jy]*(1.-b[ix,jy]**2)/6.)

	mass[p[ix,jy]+2,q[ix,jy]+2] = mass[p[ix,jy]+2,q[ix,jy]+2]-(a[ix,jy]*(1.-a[ix,jy]**2)*(a16[ix,jy]/a0)/6.)\
													*(-b[ix,jy]*(1.-b[ix,jy]**2)/6.)			

	# #Attrixbute mass to the correct departure poixnts from halovals
	# mass[nhalo:2*nhalo,:]=	mass[nhalo:2*nhalo,:]+mass[nx+nhalo:nx+2*nhalo,:]
	# mass[nx:nx+nhalo,:] = mass[nx:nx+nhalo,:]+mass[0:nhalo,:]

	# mass[:,nhalo:2*nhalo]=	mass[:,nhalo:2*nhalo]+mass[:,nx+nhalo:nx+2*nhalo]
	# mass[:,nx:nx+nhalo] = mass[:,nx:nx+nhalo]+mass[:,0:nhalo]
	
	
	#centers
	mass[nhalo:2*nhalo,nhalo:nx+nhalo] = mass[nhalo:2*nhalo,nhalo:nx+nhalo]+mass[nx+nhalo:nx+2*nhalo,nhalo:nx+nhalo]
	mass[nx:nx+nhalo,nhalo:nx+nhalo] = mass[nx:nx+nhalo,nhalo:nx+nhalo]+mass[0:nhalo,nhalo:nx+nhalo]

	mass[nhalo:nx+nhalo,nhalo:2*nhalo] = mass[nhalo:nx+nhalo,nhalo:2*nhalo]+mass[nhalo:nx+nhalo,nx+nhalo:nx+2*nhalo]
	mass[nhalo:nx+nhalo,nx:nx+nhalo] = mass[nhalo:nx+nhalo,nx:nx+nhalo]+mass[nhalo:nx+nhalo,0:nhalo]

	#corners
	mass[nx:nx+nhalo,nx:nx+nhalo] = mass[nx:nx+nhalo,nx:nx+nhalo] + mass[0:nhalo,0:nhalo]
	mass[nx:nx+nhalo,nhalo:2*nhalo] = mass[nx:nx+nhalo,nhalo:2*nhalo] + mass[0:nhalo,nx+nhalo:nx+2*nhalo]

	mass[nhalo:2*nhalo,nx:nx+nhalo] = mass[nhalo:2*nhalo,nx:nx+nhalo] + mass[nx+nhalo:nx+2*nhalo,0:nhalo]
	mass[nhalo:2*nhalo,nhalo:2*nhalo] = mass[nhalo:2*nhalo,nhalo:2*nhalo] + mass[nx+nhalo:nx+2*nhalo,nx+nhalo:nx+2*nhalo]
	


	# #W weixght
	w[nhalo:nx+nhalo,nhalo:nx+nhalo] = 1.0/mass[nhalo:nx+nhalo,nhalo:nx+nhalo]
	haloval(w)

	wa1[ix,jy] =    -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]-1]*w[p[ix,jy]-1,q[ix,jy]-1]*(a1[ix,jy]/a0)/6. \
			   +(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]-1]*w[p[ix,jy],q[ix,jy]-1]*(a2[ix,jy]/a0)/2. \
			+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]-1]*w[p[ix,jy]+1,q[ix,jy]-1]*(a3[ix,jy]/a0)/2. \
					-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]-1]*w[p[ix,jy]+2,q[ix,jy]-1]*(a4[ix,jy]/a0)/6.

	wa2[ix,jy] =    -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]]*w[p[ix,jy]-1,q[ix,jy]  ]*(a5[ix,jy]/a0)/6. \
				+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]]*w[p[ix,jy]  ,q[ix,jy]]*(a6[ix,jy]/a0)/2. \
			+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]]*w[p[ix,jy]+1,q[ix,jy]]*(a7[ix,jy]/a0)/2. \
					-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]]*w[p[ix,jy]+2,q[ix,jy]]*(a8[ix,jy]/a0)/6.

	wa3[ix,jy] =    -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+1]*w[p[ix,jy]-1,q[ix,jy]+1]*(a9[ix,jy]/a0)/6. \
				+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+1]*w[p[ix,jy],q[ix,jy]+1]*(a10[ix,jy]/a0)/2. \
			+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+1]*w[p[ix,jy]+1,q[ix,jy]+1]*(a11[ix,jy]/a0)/2. \
					-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+1]*w[p[ix,jy]+2,q[ix,jy]+1]*(a12[ix,jy]/a0)/6.
 
	wa4[ix,jy] = 	-(a[ix,jy])*(1.0-a[ix,jy])*(2.0-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+2]*w[p[ix,jy]-1,q[ix,jy]+2]*(a13[ix,jy]/a0)/6.0 \
				+(1.0-a[ix,jy]**2)*(2.0-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+2]*w[p[ix,jy],q[ix,jy]+2]*(a14[ix,jy]/a0)/2.0 \
			+a[ix,jy]*(1.0+a[ix,jy])*(2.0-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+2]*w[p[ix,jy]+1,q[ix,jy]+2]*(a15[ix,jy]/a0)/2.0 \
					-a[ix,jy]*(1.-a[ix,jy]**2) *arr[p[ix,jy]+2,q[ix,jy]+2]*w[p[ix,jy]+2,q[ix,jy]+2]*(a16[ix,jy]/a0)/6.0
					
	arr2[ix,jy] = -(b[ix,jy])*(1.0-b[ix,jy])*(2.0-b[ix,jy])*wa1[ix,jy]/6.0 \
				+(1.0-b[ix,jy]**2)*(2.0-b[ix,jy])*wa2[ix,jy]/2.0 \
				+b[ix,jy]*(1.0+b[ix,jy])*(2.0-b[ix,jy])*wa3[ix,jy]/2.0 \
					-b[ix,jy]*(1.0-b[ix,jy]**2)*wa4[ix,jy]/6.0

	haloval(arr2)
	return arr2

def bicubich(a,b,p,q,ix,jy,arr,arrA,arrB,arrC,arrD):
	#A height
	ix+=-1
	wa1[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]-1]/6. \
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]-1]/2. \
							+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]-1]/2. \
							-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]-1]/6.
						
	wa2[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]]/6.\
							+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]]/2.\
							+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]]/2.\
							-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]]/6.

	wa3[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+1]/6.\
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+1]/2.\
						+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+1]/2.\
						-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+1]/6.

	wa4[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+2]/6.\
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+2]/2.\
						+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+2]/2.\
						-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+2]/6.

	arrA[ix+1,jy] = -(b[ix,jy])*(1.-b[ix,jy])*(2.-b[ix,jy])*wa1[ix,jy]/6.\
							+(1.-b[ix,jy]**2)*(2.-b[ix,jy])*wa2[ix,jy]/2.\
							+b[ix,jy]*(1.+b[ix,jy])*(2.-b[ix,jy])*wa3[ix,jy]/2.\
							-b[ix,jy]*(1.-b[ix,jy]**2)*wa4[ix,jy]/6.			
	ix+=1
	
	
	#B height
	ix+=1
	wa1[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]-1]/6. \
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]-1]/2. \
							+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]-1]/2. \
							-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]-1]/6.
						
	wa2[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]]/6.\
							+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]]/2.\
							+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]]/2.\
							-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]]/6.

	wa3[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+1]/6.\
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+1]/2.\
						+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+1]/2.\
						-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+1]/6.

	wa4[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+2]/6.\
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+2]/2.\
						+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+2]/2.\
						-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+2]/6.

	arrB[ix-1,jy] = -(b[ix,jy])*(1.-b[ix,jy])*(2.-b[ix,jy])*wa1[ix,jy]/6.\
							+(1.-b[ix,jy]**2)*(2.-b[ix,jy])*wa2[ix,jy]/2.\
							+b[ix,jy]*(1.+b[ix,jy])*(2.-b[ix,jy])*wa3[ix,jy]/2.\
							-b[ix,jy]*(1.-b[ix,jy]**2)*wa4[ix,jy]/6.			
	ix+=-1
	


	#C height
	jy+=-1
	wa1[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]-1]/6. \
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]-1]/2. \
							+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]-1]/2. \
							-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]-1]/6.
						
	wa2[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]]/6.\
							+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]]/2.\
							+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]]/2.\
							-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]]/6.

	wa3[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+1]/6.\
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+1]/2.\
						+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+1]/2.\
						-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+1]/6.

	wa4[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+2]/6.\
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+2]/2.\
						+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+2]/2.\
						-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+2]/6.

	arrC[ix,jy+1] = -(b[ix,jy])*(1.-b[ix,jy])*(2.-b[ix,jy])*wa1[ix,jy]/6.\
							+(1.-b[ix,jy]**2)*(2.-b[ix,jy])*wa2[ix,jy]/2.\
							+b[ix,jy]*(1.+b[ix,jy])*(2.-b[ix,jy])*wa3[ix,jy]/2.\
							-b[ix,jy]*(1.-b[ix,jy]**2)*wa4[ix,jy]/6.			
	jy+=1
	
	
	
	#D height
	jy+=1
	wa1[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]-1]/6. \
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]-1]/2. \
							+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]-1]/2. \
							-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]-1]/6.
						
	wa2[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]]/6.\
							+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]]/2.\
							+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]]/2.\
							-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]]/6.

	wa3[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+1]/6.\
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+1]/2.\
						+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+1]/2.\
						-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+1]/6.

	wa4[ix,jy] = -(a[ix,jy])*(1.-a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]-1,q[ix,jy]+2]/6.\
						+(1.-a[ix,jy]**2)*(2.-a[ix,jy])*arr[p[ix,jy],q[ix,jy]+2]/2.\
						+a[ix,jy]*(1.+a[ix,jy])*(2.-a[ix,jy])*arr[p[ix,jy]+1,q[ix,jy]+2]/2.\
						-a[ix,jy]*(1.-a[ix,jy]**2)*arr[p[ix,jy]+2,q[ix,jy]+2]/6.

	arrD[ix,jy-1] = -(b[ix,jy])*(1.-b[ix,jy])*(2.-b[ix,jy])*wa1[ix,jy]/6.\
							+(1.-b[ix,jy]**2)*(2.-b[ix,jy])*wa2[ix,jy]/2.\
							+b[ix,jy]*(1.+b[ix,jy])*(2.-b[ix,jy])*wa3[ix,jy]/2.\
							-b[ix,jy]*(1.-b[ix,jy]**2)*wa4[ix,jy]/6.			
	jy+=-1
	
	haloval(arrA)
	haloval(arrB)
	haloval(arrC)
	haloval(arrD)
	return
	
nx = 128												#gridpoints x
ny = nx													#gridpoints y
nhalo = 10												#gridpoints to make grid periodic
xmin = 0.0
xmax = 1000000.0
Lx = xmax-xmin
dx = Lx/(nx-1)
ymin = xmin
ymax = xmax
Ly = ymax-ymin
dy = Ly/(ny-1)
H0 = 4000.0
nstop = 30
u0 = 10.0 #30
v0 = u0
g = 9.82
cg = sqrt(g*H0)
tau = (1.0/sqrt(2.0))*dx/(sqrt(u0**2+v0**2)+cg)
Ns = 10
f = 10**(-4) #coriolis parameter
a0 = 0.25*2.0*dx*2.0*dy

#Tstop = 30.0*Lx/cg
#tau1 = 2.0*dx/cg
#mult = 1.5
#dt1 = mult*Ns*tau1
#nstop=floor(Tstop/dt1)+1
#dt = Tstop/nstop
#tau = dt/Ns


#Filter coefficients
w1 = 1.0/10.0
w0 = 1.0 - 4.0*w1
#offcentering variable
gamma = 1.05
#Iteration of trajectories, guessing departure point
niter = 3

x = linspace(0,Lx,nx)

datatype=float64


h = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype) #height field
d = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype)  #depth, h-hs
hstilde = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype) #short timestep bottom topy
hs = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype) #bottom topography
hA = zeros((nx+2*nhalo,nx+2*nhalo))
hB = zeros((nx+2*nhalo,nx+2*nhalo))
hC = zeros((nx+2*nhalo,nx+2*nhalo))
hD = zeros((nx+2*nhalo,nx+2*nhalo))
Tx = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype) #downstream trajectory
Ty = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype) #downstream trajectory
Cupx = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype) #upstream trajectory for final advection
Cupy = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype) #upstream trajectory for final advection
u = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype) #horizontal velocity
u[0,:,:] = array(u[0,:,:]+u0)
u[1,:,:] = u[0,:,:]
v = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype) #horizontal velocity
v[0,:,:] = array(v[0,:,:]+u0)
v[1,:,:] = v[0,:,:]

a = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
p = zeros((nx+2*nhalo,nx+2*nhalo),dtype=int) #int or int8
b = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
q = zeros((nx+2*nhalo,nx+2*nhalo),dtype=int)
ii = zeros((nx+2*nhalo,nx+2*nhalo),dtype=int)
jj = zeros((nx+2*nhalo,nx+2*nhalo),dtype=int)

#xdep = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
#ydep = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)

work1 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
work2 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
work2x = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)# just use work1 and work 2, efficiently..
work2y = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
area1 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype) #Can perhaps again here, use work1 and work2
area2 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
arr = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
arr2 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
mass = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype) #maybe just use work1
w = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype) #maybe just use work2
#w = zeros((nx+2*nhalo,nx+2*nhalo))
wa1 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
wa2 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
wa3 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
wa4 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
#areas used 
a1 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a2 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a3 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a4 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a5 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a6 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a7 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a8 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a9 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a10 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a11 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a12 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a13 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a14 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a15 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)
a16 = zeros((nx+2*nhalo,nx+2*nhalo),dtype=dtype)


#Initialize the grid point index vector, ii,jj
for i in range(nhalo,nx+nhalo): 
  ii[i,:] = int(i)#int(i-nhalo)
  jj[:,i] = int(i)#int(i-nhalo)
haloval(ii)
haloval(jj)
#vectorization indices
ix = ii[nhalo:nx+nhalo,nhalo:nx+nhalo]
jy = jj[nhalo:nx+nhalo,nhalo:nx+nhalo]

#bottom topography
#for i in range(nhalo,nx+nhalo):
#	for j in range(nhalo,nx+nhalo):
#	r = abs(5.0*2.0*pi*(x[i-nhalo]-0.5*(xmin+xmax))/(xmax-xmin))
#	r = min(pi,r)
#	hs[i] = 0.1*H0*(1.0+cos(r))
#		hs[i,j] = 0.0001*H0
#hs = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
hs[ix,jy] = 0.1*H0*exp(-(ix-nx/2)**2/(nx/3)-(jy-ny/2)**2/(ny/3)) #Gauss mountain
haloval(hs[:,:])


#Initial h field
#for i in range(nhalo,nx+nhalo):
#	for j in range(nhalo,nx+nhalo):
#		h[0,i,j] = H0+1.0*exp(-(i-nx/2)**2/(nx/3)-(j-ny/2)**2/(ny/3)) #exponential wave
h[0,ix,jy] = H0+1.0*exp(-(ix-nx/2)**2/(nx/3)-(jy-ny/2)**2/(ny/3)) #Gauss wave

haloval(h[0,:,:])
h[1,:,:]=h[0,:,:]
h[2,:,:]=h[0,:,:]
h[3,:,:]=h[0,:,:]

#Initial d field
d[0,:,:] = h[0,:,:]-hs
haloval(d[0,:,:])
#d[1,:,:] = h[0,:,:]-hs[:,:]
#d[2,:,:] = h[0,:,:]-hs[:,:]
#d[3,:,:] = h[0,:,:]-hs[:,:]


n1 = 0
n2 = 1
ns1 = 2
ns2 = 3


#saving bottom topo and initial state
for i in range(nhalo,nx+nhalo):
	for j in range(nhalo,ny+nhalo):
		f3.write('{0:30.15f}'.format(hs[i,j]))
	f3.write("\n")
#savetxt("bottomtopo2d.txt",hs[nhalo:nx+nhalo,nhalo:nx+nhalo])	
#save("bottomtopo2d.npy",hs[nhalo:nx+nhalo,nhalo:nx+nhalo])	

for i in range(nhalo,nx+nhalo):
	for j in range(nhalo,ny+nhalo):
		f2.write('{0:20.5f}'.format(h[0,i,j]))
	f2.write("\n")
#savetxt("data2d.txt",h[0,nhalo:nx+nhalo,nhalo:nx+nhalo])
#save("data2d.npy",h[0,nhalo:nx+nhalo,nhalo:nx+nhalo])
#savetxt(f2,h[0,nhalo:nx+nhalo,nhalo:nx+nhalo])
#savetxt(f2, h[0,nhalo:nx+nhalo,nhalo:nx+nhalo])

#http://stackoverflow.com/questions/27786868/python3-numpy-appending-to-a-file-using-numpy-savetxt
#savetxt(f2,h[0,nhalo:nx+nhalo,nhalo:nx+nhalo])

#total mass
summass = d[0,nhalo:nx+nhalo,nhalo:nx+nhalo].sum()
#Saving total initial mass
f4.write('{0:20.30f}'.format(summass))
f4.write("\n")

print('n={0:3.0f} tau={1:4.5f} hmax={2:4.2f} umax={3:4.2f} umin={4:4.2f} vmax={5:4.2f} vmin={6:4.2f} mass={7:10.2f}'\
		.format(0, tau, h[n2,:,:].max(), u[n2,:,:].max(), u[n2,:,:].min(), v[n2,:,:].max(), v[n2,:,:].min(), summass))

for n in range(int(nstop)): #int(nstop)
	#We need to set T = 0 for both short timestep pointers, or well, AT LEAST, for T[2,:] = 0
	#Tx[2,:,:] = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
	#Tx[3,:,:] = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
	Tx = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype)
	Ty = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype)
	#Ty[2,:,:] = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
	#Ty[3,:,:] = zeros((nx+2*nhalo,nx+2*nhalo),dtype=datatype)
	Cupx = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype)
	Cupy = zeros((4,nx+2*nhalo,nx+2*nhalo),dtype=datatype)
	
	#short timestep pointers
	ns1 = 2
	ns2 = 3
	
	#We update with the newly advected velocities, affected by iteration of trajectories
	u[2,:,:] = u[n1,:,:]
	v[2,:,:] = v[n1,:,:]
	h[2,:,:] = h[n1,:,:]
	d[2,:,:] = d[n1,:,:]
	#It's already been halovalled, so no need to call again
	#hstilde[0,:,:] = hs
	#hstilde[1,:,:] = hs
	hstilde[2,:,:] = hs
	#hstilde[3,:,:] = hs


	#Short timeloop, forward-backward integration
	for ns in range(Ns):#Ns
		#updating short timestep tau, based on u,v and h.
		maxvalu = absolute(u[ns1,nhalo:nx+nhalo,nhalo:nx+nhalo]).max()
		maxvalv = absolute(v[ns1,nhalo:nx+nhalo,nhalo:nx+nhalo]).max()
		maxvalcg = sqrt(g*h[ns1,nhalo:nx+nhalo,nhalo:nx+nhalo].max())
		tau = (1.0/sqrt(2.0))*dx/(sqrt(maxvalu**2+maxvalv**2)+cg)					
		

		work1 = gamma*h[ns1,:,:]+(1.0-gamma)*h[ns2,:,:]
		haloval(work1)
		
		#Arrival point, no minus sign.
		#xdep[:,:] = Tx[ns1,:,:]/dx 
		p[:,:] = floor(Tx[ns1,:,:]/dx)
		a[:,:] = Tx[ns1,:,:]/dx-p[:,:]
		p[:,:] += ii[:,:]
	
		#ydep[:,:] = Ty[ns1,:,:]/dy
		q[:,:] = floor(Ty[ns1,:,:]/dy)
		b[:,:] = Ty[ns1,:,:]/dy-q[:,:]
		q[:,:] += jj[:,:]
	
		haloval(a)
		haloval(p)
		haloval(b)
		haloval(q)

		haloval(h[ns1,:,:])
		
		#bilinearh(a,b,p,q,ix,jy,work1,hA,hB,hC,hD)
		bicubich(a,b,p,q,ix,jy,work1,hA,hB,hC,hD)

		#modified forward velocity
		work1[nhalo:nx+nhalo,:] = -tau*g*(hB[nhalo:nx+nhalo,:]-hA[nhalo:nx+nhalo,:])\
									/(2.0*dx + Tx[ns1,nhalo+1:nx+nhalo+1,:] - Tx[ns1,nhalo-1:nx+nhalo-1,:])

		work2[:,nhalo:ny+nhalo] = -tau*g*(hD[:,nhalo:ny+nhalo]-hC[:,nhalo:ny+nhalo])\
									/(2.0*dy + Ty[ns1,:,nhalo+1:ny+nhalo+1] - Ty[ns1,:,nhalo-1:ny+nhalo-1])

		#Normal forward velocity					
		#work2x[nhalo:nx+nhalo,:] = -tau*g*(work1[nhalo+1:nx+nhalo+1,:]-work1[nhalo-1:nx+nhalo-1,:])\
		#							 /(2.0*dx + Tx[ns1,nhalo+1:nx+nhalo+1,:] - Tx[ns1,nhalo-1:nx+nhalo-1,:])

		#work2y[:,nhalo:ny+nhalo] = -tau*g*(work1[:,nhalo+1:ny+nhalo+1]-work1[:,nhalo-1:ny+nhalo-1])\
		#							 /(2.0*dy + Ty[ns1,:,nhalo+1:ny+nhalo+1] - Ty[ns1,:,nhalo-1:ny+nhalo-1])
	
		#forward forecast velocity
		u[ns2,:,:] = u[ns1,:,:]+work1[:,:]+tau*0.0001*v[ns1,:,:]
		v[ns2,:,:] = v[ns1,:,:]+work2[:,:]-tau*0.0001*u[ns1,:,:]
		haloval(u[ns2,:,:])
		haloval(v[ns2,:,:])
		#trajectory
		Tx[ns2,:,:] = Tx[ns1,:,:]+tau*u[ns2,:,:]
		Ty[ns2,:,:] = Ty[ns1,:,:]+tau*v[ns2,:,:]
		haloval(Tx[ns2,:,:])
		haloval(Ty[ns2,:,:])		



		# Backwards forecast of mass, by depth:	
		#These are quasi-lagrangian areas, deformed by trajectories. Approximate areas.
		area1[nhalo:nx+nhalo,nhalo:nx+nhalo] = sqrt((2.0*dx + Tx[ns1,nhalo+1:nx+nhalo+1,nhalo:nx+nhalo] - Tx[ns1,nhalo-1:nx+nhalo-1,nhalo:nx+nhalo])**2\
					 +(Ty[ns1,nhalo+1:nx+nhalo+1,nhalo:nx+nhalo]-Ty[ns1,nhalo-1:nx+nhalo-1,nhalo:nx+nhalo])**2)\
					 *sqrt((2.0*dx + Ty[ns1,nhalo:nx+nhalo,nhalo+1:nx+nhalo+1] - Ty[ns1,nhalo:nx+nhalo,nhalo-1:nx+nhalo-1])**2\
					 +(Tx[ns1,nhalo:nx+nhalo,nhalo+1:nx+nhalo+1]-Tx[ns1,nhalo:nx+nhalo,nhalo-1:nx+nhalo-1])**2)
	
		area2[nhalo:nx+nhalo,nhalo:nx+nhalo] =  sqrt((2.0*dx + Tx[ns2,nhalo+1:nx+nhalo+1,nhalo:nx+nhalo] - Tx[ns2,nhalo-1:nx+nhalo-1,nhalo:nx+nhalo])**2\
					 +(Ty[ns2,nhalo+1:nx+nhalo+1,nhalo:nx+nhalo]-Ty[ns2,nhalo-1:nx+nhalo-1,nhalo:nx+nhalo])**2)\
					 *sqrt((2.0*dx + Ty[ns2,nhalo:nx+nhalo,nhalo+1:nx+nhalo+1] - Ty[ns2,nhalo:nx+nhalo,nhalo-1:nx+nhalo-1])**2\
					 +(Tx[ns2,nhalo:nx+nhalo,nhalo+1:nx+nhalo+1]-Tx[ns2,nhalo:nx+nhalo,nhalo-1:nx+nhalo-1])**2)
			
		d[ns2,nhalo:nx+nhalo,nhalo:nx+nhalo] = d[ns1,nhalo:nx+nhalo,nhalo:nx+nhalo]\
											*area1[nhalo:nx+nhalo,nhalo:nx+nhalo]/area2[nhalo:nx+nhalo,nhalo:nx+nhalo]

		haloval(d[ns2,:,:])
		
		#Downstream hs arrival
		#Arrival point, no minus sign.
		#xdep[:,:] = Tx[ns2,:,:]/dx
		p[:,:] = floor(Tx[ns2,:,:]/dx)
		a[:,:] = Tx[ns2,:,:]/dx-p[:,:]
		p[:,:] += ii[:,:]
	
		#ydep[:,:] = Ty[ns2,:,:]/dy
		q[:,:] = floor(Ty[ns2,:,:]/dy)
		b[:,:] = Ty[ns2,:,:]/dy-q[:,:]
		q[:,:] += jj[:,:] #do in-place?
	
		haloval(a)
		haloval(p)
		haloval(b)
		haloval(q)
		
		#Interpolation of hs at arrival point
		bicubicvec(a,b,p,q,hs,hstilde[ns2,:,:])
		
		#total height
		h[ns2,:,:] = d[ns2,:,:]+hstilde[ns2,:,:]
		haloval(h[ns2,:,:])
	
		#Filtering waves
		work2[:,:] = h[ns2,:,:] - h[ns1,:,:]
		for i in range(2):
			work1[nhalo:nx+nhalo,nhalo:nx+nhalo] = w1*work2[nhalo-1:nx+nhalo-1,nhalo:nx+nhalo] + w0*work2[nhalo:nx+nhalo,nhalo:nx+nhalo]\
								 + w1*work2[nhalo+1:nx+nhalo+1,nhalo:nx+nhalo]\
								 +w1*work2[nhalo:nx+nhalo,nhalo-1:nx+nhalo-1] + w1*work2[nhalo:nx+nhalo,nhalo+1:nx+nhalo+1]
			haloval(work1)
		
			work2[nhalo:nx+nhalo,nhalo:nx+nhalo] = w1*work1[nhalo-1:nx+nhalo-1,nhalo:nx+nhalo] + w0*work1[nhalo:nx+nhalo,nhalo:nx+nhalo]\
								 +w1*work1[nhalo+1:nx+nhalo+1,nhalo:nx+nhalo]\
                                +w1*work1[nhalo:nx+nhalo,nhalo-1:nx+nhalo-1] + w1*work1[nhalo:nx+nhalo,nhalo+1:nx+nhalo+1]
			haloval(work2)

			
		#update filtered height	
		h[ns2,:,:] = h[ns1,:,:]+work2[:,:]
		
		#after filtering heights
		haloval(h[ns2,:,:])
		
		#update depth based on filtered heights
		d[ns2,:,:] = h[ns2,:,:] - hstilde[ns2,:,:]
		haloval(d[ns2,:,:])
		
		#changing pointers
		#nst = ns1
		#ns1 = ns2
		#ns2 = nst
		ns1,ns2 = ns2,ns1
	
	#T[2,:] = T[ns1,:]
	u[2,:,:] = u[ns1,:,:]
	v[2,:,:] = v[ns1,:,:]
	h[2,:,:] = h[ns1,:,:] 
	d[2,:,:] = d[ns1,:,:] #Maybe do not update this. Well it doesn't matter
	#Since i also have d[2,:] = d[n1,:] right before short timestep, so hm..
	
	
	#Shouldn't be necessary...
	# haloval(u[2,:,:])
	# haloval(v[2,:,:])
	# haloval(h[2,:,:])
	# haloval(d[2,:,:])
	# haloval(Tx[ns1,:,:])
	# haloval(Ty[ns1,:,:])
	
	#first guess departure point
	#xdep[:,:] = -Tx[ns1,:,:]/dx
	p[:,:] = floor(-Tx[ns1,:,:]/dx)
	a[:,:] = -Tx[ns1,:,:]/dx-p[:,:]
	p[:,:] += ii[:,:]
	
	#ydep[:,:] = -Ty[ns1,:,:]/dy
	q[:,:] = floor(-Ty[ns1,:,:]/dy)
	b[:,:] = -Ty[ns1,:,:]/dy-q[:,:]
	q[:,:] += jj[:,:]
	
	haloval(a)
	haloval(p)
	haloval(b)
	haloval(q)
	
	#first guess trajectory
	bicubicvec(a,b,p,q,Tx[ns1,:,:],Cupx[ns1,:,:])
	bicubicvec(a,b,p,q,Ty[ns1,:,:],Cupy[ns1,:,:])

	
	for it in range(niter):
		#xdep[:,:] = -Cupx[ns1,:,:]/dx
		p[:,:] = floor(-Cupx[ns1,:,:]/dx)
		a[:,:] = -Cupx[ns1,:,:]/dx-p[:,:]
		p[:,:] += ii[:,:]
	
		#ydep[:,:] = -Cupy[ns1,:,:]/dy
		q[:,:] = floor(-Cupy[ns1,:,:]/dy)
		b[:,:] = -Cupy[ns1,:,:]/dy-q[:,:]
		q[:,:] += jj[:,:]
	
		haloval(a)
		haloval(p)
		haloval(b)
		haloval(q)
		#cubicint(a,p,T[ns1,:],Cup[ns1,:])
		bicubicvec(a,b,p,q,Tx[ns1,:,:],Cupx[ns1,:,:])
		bicubicvec(a,b,p,q,Ty[ns1,:,:],Cupy[ns1,:,:])
	
	#Upstream trajectory to a T landing in grid proint, for final advection of fields
	#xdep[:,:] = -Cupx[ns1,:,:]/dx
	p[:,:] = floor(-Cupx[ns1,:,:]/dx)
	a[:,:] = -Cupx[ns1,:,:]/dx-p[:,:]
	p[:,:] += ii[:,:]
	
	#ydep[:,:] = -Cupy[ns1,:,:]/dy
	q[:,:] = floor(-Cupy[ns1,:,:]/dy)
	b[:,:] = -Cupy[ns1,:,:]/dy-q[:,:]
	q[:,:] += jj[:,:]	
	
	haloval(a)
	haloval(p)
	haloval(b)
	haloval(q)
	#Advection of velocity
	#cubicint(a[:],p[:],u[2,:],u[n2,:])
	bicubicvec(a,b,p,q,u[2,:,:],u[n2,:,:])
	bicubicvec(a,b,p,q,v[2,:,:],v[n2,:,:])
	
	#Advection of height field
	#bilin(a,b,p,q,h[2,:,:],h[n2,:,:])
	#bicubicvec(a,b,p,q,h[2,:,:],h[n2,:,:])
	#update depth
	#d[n2,:,:] = h[n2,:,:]-hs
	#haloval(d[n2,:,:])
	
	
	#Advection of d field
	#bicubic(a,b,p,q,d[2,:,:],d[n2,:,:])
	bicubicvec(a,b,p,q,d[2,:,:],d[n2,:,:])

	#update height

	#LMCSL, from n1 to n2.
	#bicubicLMCSL(a,b,p,q,Cupx[ns1,:,:],Cupy[ns1,:,:],d[n1,:,:],d[n2,:,:],a0,w)
	#bicubicLMCSLvec(a,b,p,q,Cupx[ns1,:,:],Cupy[ns1,:,:],d[n1,:,:],d[n2,:,:])
	#update height
	h[n2,:,:] = d[n2,:,:]+hs


	#Filtering waves
	work2[:,:] = h[n2,:,:] - h[n1,:,:]
	for i in range(2):
		work1[nhalo:nx+nhalo,nhalo:nx+nhalo] = w1*work2[nhalo-1:nx+nhalo-1,nhalo:nx+nhalo] + w0*work2[nhalo:nx+nhalo,nhalo:nx+nhalo]\
							 + w1*work2[nhalo+1:nx+nhalo+1,nhalo:nx+nhalo]\
							 +w1*work2[nhalo:nx+nhalo,nhalo-1:nx+nhalo-1] + w1*work2[nhalo:nx+nhalo,nhalo+1:nx+nhalo+1]
		haloval(work1)
	
		work2[nhalo:nx+nhalo,nhalo:nx+nhalo] = w1*work1[nhalo-1:nx+nhalo-1,nhalo:nx+nhalo] + w0*work1[nhalo:nx+nhalo,nhalo:nx+nhalo]\
							 +w1*work1[nhalo+1:nx+nhalo+1,nhalo:nx+nhalo]\
							+w1*work1[nhalo:nx+nhalo,nhalo-1:nx+nhalo-1] + w1*work1[nhalo:nx+nhalo,nhalo+1:nx+nhalo+1]
		haloval(work2)

			
	#update filtered height	
	h[n2,:,:] = h[n1,:,:]+work2[:,:]
	
	#after filtering heights.
	haloval(h[n2,:,:])

	d[n2,:,:] = h[n2,:,:]-hs[:,:]
	haloval(d[n2,:,:])

	
	summass = d[n2,nhalo:nx+nhalo,nhalo:nx+nhalo].sum()
	#print('n=%s tau = %4.5f hmax=%4.2f umax=%4.2f umin=%4.2f mass=%10.2f'\
	#		% (n+1, tau, h[n2,:,:].max(), u[n2,:,:].max(), u[n2,:,:].min(), summass))
			
	print('n={0:3.0f} tau={1:4.5f} hmax={2:4.2f} umax={3:4.2f} umin={4:4.2f} vmax={5:4.2f} vmin={6:4.2f} mass={7:10.2f}'\
		.format(n+1, tau, h[n2,:,:].max(), u[n2,:,:].max(), u[n2,:,:].min(), v[n2,:,:].max(), v[n2,:,:].min(), summass))
	#savetxt(f,h[n,nhalo:nx+nhalo],fmt='%4.4f', newline="\r\n",delimiter=" ")		
	#savetxt(f,h[n,nhalo:nx+nhalo],fmt='%4.4f', newline="\n",delimiter=" ")
	
	
	#Saving total initial mass
	f4.write('{0:20.30f}'.format(summass))
	f4.write("\n")
	
	#Saving height data
	for i in range(nhalo,nx+nhalo):
		for j in range(nhalo,ny+nhalo):
			f2.write('{0:20.5f}'.format(h[n2,i,j]))
		f2.write("\n")
	
	#savetxt("data2d.txt",h[n2,nhalo:nx+nhalo,nhalo:nx+nhalo])
	#save("data2d.npy",h[n2,nhalo:nx+nhalo,nhalo:nx+nhalo],'a') #'a' or no 'a'
	#savetxt(f2,h[n2,nhalo:nx+nhalo,nhalo:nx+nhalo])#virker heller ikke
	
	#De skriver måske tilføj "\n"
	#f2.write("\n")
	#savetxt(f2, h[n2,nhalo:nx+nhalo,nhalo:nx+nhalo])
	
	#http://stackoverflow.com/questions/27786868/python3-numpy-appending-to-a-file-using-numpy-savetxt
	#savetxt(f2,h[n2,nhalo:nx+nhalo,nhalo:nx+nhalo])
	#f2.write("\n") #Can't write this, can't write str... ?
	#nt = n1
	#n1 = n2
	#n2 = nt
	n1, n2 = n2,n1
	
#print(a)
#print(xdep[:])
#print(p)
print('----------------PARAMETERS-------------------------')
print('end tau = {0:4.2f}'.format(tau))
print('runtime = {0:4.2f}'.format(time.process_time()-start_time))


f1.write('{0:4.1f} {1:4.1f} {2:4.1f} {3:4.1f} {4:4.1f} {5:4.1f} {6:4.1f}'\
		.format(nstop+1,dx,nx,Lx,ny,Ly,H0))
#savetxt(f1,[nt,dx,nx,Lx])

f1.close()
f2.close()
f3.close()
f4.close()

