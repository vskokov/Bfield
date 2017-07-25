import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from scipy import integrate
from scipy import stats
import random 

import sys

b = float(sys.argv[1])

import pylab as py

R_Cu=4.2 
a_Cu=0.596
Z_Cu=29

R_Au=6.38
a_Au=0.535
Z_Au=79

R=R_Au
a=a_Au
Z=Z_Au

Rcut=0.3

x0=0.0
y0=0.0
z0=0.0

sqrts=200.0

Nevent=200000

def density(r):
	return 1.0/(1.0+np.exp((r-R)/a))



def fermi(length):
	x=[]
	y=[]
	z=[]
	while (len(x)<length):
		xx = 2.0*R*random.uniform(-1.0, 1.0)
		yy = 2.0*R*random.uniform(-1.0, 1.0)
		zz = 2.0*R*random.uniform(-1.0, 1.0)
		r=np.sqrt(xx**2+yy**2+zz**2)
		#r=2.0*R_Cu*random.uniform(0.0, 1.0)
		rho_sample=density(0.0)*random.uniform(0.0, 1.0)
		rho=density(r)
		if rho_sample<rho:
			x.append(xx)
			y.append(yy)
			z.append(zz)
	return (np.array(x),np.array(y),np.array(z))

alpha=1.0/137.0

def velocity(sqrts):
	return(np.sqrt(1.0-4.0/sqrts**2)	)

vz=velocity(sqrts)

def magnetic_field(x,y,z,sign):
	Rx = x-x0
	Ry = y-y0
	Rz = z-z0
	Rperp2=(Rx**2+Ry**2)
	R=np.sqrt(Rperp2+Rz**2)
	prefactor=sign*alpha*vz*(vz**2-1.0)
	denom=np.sqrt(1.0-vz**2*Rperp2/R**2)**3
	Bx=prefactor/(R**3*denom)*Ry
	By=-prefactor/(R**3*denom)*Rx
	if (R<Rcut):
		Bx=0.0
		By=0.0
	return (Bx,By) 

def electric_field(x,y,z,sign):
	Rx = x-x0
	Ry = y-y0
	Rz = z-z0
	Rperp2=(Rx**2+Ry**2)
	R=np.sqrt(Rperp2+Rz**2)
	prefactor=alpha*(vz**2-1.0)
	denom=np.sqrt(1.0-vz**2*Rperp2/R**2)**3
	Ex=prefactor/(R**3*denom)*Rx
	Ey=prefactor/(R**3*denom)*Ry
	if (R<Rcut):
		Ex=0.0
		Ey=0.0
	return (Ex,Ey) 



for event in range(Nevent):
	distrT=fermi(Z)
	xT = distrT[0]+b/2.0
	yT = distrT[1]
	zT = distrT[2]/(sqrts/2.0)

	distrP=fermi(Z)
	xP = distrP[0]-b/2.0
	yP = distrP[1]
	zP = distrP[2]/(sqrts/2.0)

	Bx=0.0
	By=0.0
	Ex=0.0
	Ey=0.0
	for i in range(Z):
		Bfield=magnetic_field(xT[i],yT[i],zT[i],+1.0)
		Bx=Bx+Bfield[0]
		By=By+Bfield[1]

		Bfield=magnetic_field(xP[i],yP[i],zP[i],-1.0)
		Bx=Bx+Bfield[0]
		By=By+Bfield[1]

		Efield=electric_field(xT[i],yT[i],zT[i],+1.0)
		Ex=Ex+Efield[0]
		Ey=Ey+Efield[1]

		Efield=electric_field(xP[i],yP[i],zP[i],-1.0)
		Ex=Ex+Efield[0]
		Ey=Ey+Efield[1]

	print event, Bx*(197.327/140.0)**2, By*(197.327/140.0)**2, Ex*(197.327/140.0)**2, Ey*(197.327/140.0)**2
