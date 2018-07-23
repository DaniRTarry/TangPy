"""
@author: Tarry

type: FUNCTION

description: Function that contains all the physical parametrization of the Tang Model

input: x1,x2,time,parameters
output: psi,ug,vg,ua,va,divergence 

"""
import numpy as np
#import math as mt

### Define the function:
def tm_equations(x1,x2,x3,t,parameters):

	###-------------------
	### Parametrization
	###-------------------

	if parameters == 0:
		KK = 427e3		# Wavelength in the x1 direction
		LL = KK			# Distance Between nodal surfaces
		N1 = 6e-3		# Brunt-Vaisala frequency of layer 1
		N2 = N1/3		# Brunt-Vaisala frequency of layer 2
		U0 = 0.6		# Basic zonal current
		h1 = 1e3		# Depth of Motion Layer (layer 1)
		h2 = 1e3		# Depth of Quiescent Layer (layer 2)
	elif parameters == 1:
		KK = 250e3		# Wavelength in the x1 direction
		LL = KK/2		# Distance Between nodal surfaces
		N1 = 5.8e-3		# Brunt-Vaisala frequency of layer 1
		N2 = 2.3e-6		# Brunt-Vaisala frequency of layer 2
		U0 = 0.3		# Basic zonal current
		h1 = 500.		# Depth of Motion Layer (layer 1)
		h2 = 2000.		# Depth of Quiescent Layer (layer 2)

	# Defined by us in order to obtained the desired magnitudes
	A = 5000.

	# f and Beta obtention
	OMEGA = 7.3e-5
	Rt = 6371e3
	LAT = 38.
	f = 2*OMEGA*np.sin(np.radians(LAT))
	BETA = 2*OMEGA*np.cos(np.radians(LAT))/Rt

	# Wave numbers:
	k = 2.0*np.pi/(KK)
	l = np.pi/(LL) ### ON ANANDA'S NOTES THIS ONE IS MULTIPLIED BY 2
	MU = np.sqrt(l**2+k**2)
	
	# Coefficient between stabilities, stability wave cutoff at R=3
	R = N1/N2
	
	k1 = MU*(h1*N1)/f
	k2 = MU*(h2*N2)/f
	XI = k1*x3/h1
	
	KKK = np.tanh(k2)/np.tanh(k1)
	
	a = 2*((k1/2)-np.tanh(k1/2))/(k1-np.tanh(k1))
	b = 2*((k1/2)-1/np.tanh(k1/2))/(k1-np.tanh(k1))

	# c is the complez phase speed:
	# c = cr + j*ci
	cr = U0/2*(1- R*KKK*np.tanh(k1)/(k1*(1+R*KKK)))
	ci = U0*(k1-np.tanh(k1))*np.sqrt(-(R*KKK + a)*(R*KKK + b))/(2*k1*(1 + R*KKK))

	# Cos variable
	X = x1-cr*t
	# Exponential variable
	NU = k*ci
	T = np.log(2)/NU/3600/24

	print('cr='+str(cr)+' m/s  ||  ci='+str(ci)+' m/s  ||  T='+str(T))

	# Real and complex amplitudes
	PSIr = A*(np.cosh(XI)+k1*((cr/U0)*(1+R*KKK)-1)*np.sinh(XI)/(k1/np.tanh(k1)-1))
	PSIi = A*k1*(ci/U0)*(1+R*KKK)*np.sinh(XI)/(k1/np.tanh(k1)-1)
	# Amplitude of the wave
	PSIabs = np.sqrt(PSIr*PSIr + PSIi*PSIi)
	THETA = np.arctan(PSIi/PSIr)

	# For the vertical movement
	Br = k1*((cr/U0)*(1+R*KKK)-1)/(k1/np.tanh(k1)-1)
	Bi = k1*(ci/U0)*(1+R*KKK)/(k1/np.tanh(k1)-1)

	Wr1 = (k1*ci/U0+Bi)*np.sinh(XI)
	Wr2 = (XI*Bi-k1/U0*(cr*Bi+ci*Br))*np.cosh(XI)
	Wr = -f*U0*k*A/N1/N1/h1*(Wr1-Wr2)

	Wi1 = (XI-k1*cr/U0-Br)*np.sinh(XI)
	Wi2 = (XI*Br-1-k1/U0*(cr*Br-ci*Bi))*np.cosh(XI)
	Wi = -f*U0*k*A/N1/N1/h1*(Wi1+Wi2)

	Wabs = np.sqrt(Wr*Wr+Wi*Wi)
	SIGMA = np.arctan(-Wi/Wr)





	###-------------------
	### Model Equations
	###-------------------

	# To simplify the model equations we will calculate beforehand
	cosX = np.cos(k*X+THETA)
	sinX = np.sin(k*X+THETA)
	cosy = np.cos(l*x2)
	siny = np.sin(l*x2)
	expt = np.exp(NU*t)

	# The Stream Function
	psi = PSIabs*cosX*expt*siny

	# The geostrophic velocities
	ug = -expt*l*PSIabs*cosy*cosX
	vg = -expt*k*PSIabs*siny*sinX
	
	
	# The ageostrophic velocities
	ua1 = l*cosy*(-x2*BETA*cosX+expt*k*k*PSIabs*siny)
	ua2 = k*siny*(cr*k*cosX-NU*sinX)
	ua = -1/f*expt*PSIabs*(ua1+ua2)
	
	# va1 = -2*k*x2*BETA*siny*(-1)*sinX
	# va2 = 2*l*cosy*(-NU*cosX+cr*k*(-1)*sinX)
	# va3 = expt*k*l*l*PSIabs*np.sin(2*k*X+2*THETA)
	# va = 1/2/f*expt*PSIabs*(va1+va2+va3)
	va = (1./(2.*f))*np.exp(t*NU)*PSIabs*(-2*k*x2*BETA*np.sin(l*x2)*np.sin(cr*k*t-k*x1-THETA)+2*l*np.cos(l*x2)*(-NU*np.cos(cr*k*t-k*x1-THETA)+cr*k*np.sin(cr*k*t-k*x1-THETA))+np.exp(t*NU)*k*l*l*PSIabs*np.sin(2*cr*k*t-2*(k*x1+THETA)))


	# The horizontal divergence
	div1 = (k*k+l*l)*NU*cosX*siny
	div2 = k*(cr*(k*k+l*l)+BETA)*sinX
	divergence = 1/f*expt*PSIabs*siny*(div1+div2)


	# Vorticity
	vort1 = expt*k*k*l*l*PSIabs*(np.cos(2*l*x2)-np.cos(2*k*X+2*THETA))
	vort2 = cosX*(-1*l*BETA*cosy+(k*k+l*l)*(-f+x2*BETA)*siny)
	vorticity = 1/f*expt*PSIabs*(vort1+vort2)


	# Vertical velocity
	w = Wabs*np.cos(k*X-SIGMA)*expt*siny


	return(psi,ug,vg,ua,va,w,divergence,vorticity)



