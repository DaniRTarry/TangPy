# -*- coding: utf-8 -*-
"""
@author: Tarry

type: PROGRAM

description: Main program that generates the fields from tang's model

"""
print('Running file: tang_model.py')

### Import libraries
#import math as mt
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import os

### Import functions
from tang_model_equations import tm_equations
from tang_model_plot import tm_cfplot, tm_vec_cfplot, tm_vecplot
from gifcreator import tm_gif






###------------------------------------------------
### Define options and parameters of the program 	
###------------------------------------------------

### Options
plot_psi 	 	= False
plot_dh  	 	= True
plot_vec_Ugt 	= False
plot_vec_Uat 	= False
plot_vec_Ut  	= False
plot_Ugt_dh  	= True
plot_Uat_dh  	= True
plot_Ut_dh   	= False
plot_div_dh	 	= True
plot_Uat_div 	= False
plot_Uat_div_dh = True
plot_div_w		= False
plot_w_dh		= True
plot_vort_dh 	= True

plot_gif 		= True

# Print options
#print('plot_psi = '+str(plot_psi))
#print('plot_dh_adj = '+str(plot_dh_adj))
#print('plot_vel = '+str(plot_vel))
#print('plot_gif = '+str(plot_gif))


### Parameters:
# We need to define 2 cases:
# 0 for Tang Model parameters
# 1 for Algerian Basin parameters
parameters = 1
if parameters == 0:
	zone = 'Tang_paper'
elif parameters == 1:
	zone = 'Algerian_Basin'
print('Model Parameters set for: '+zone)

# Set Image directory
img_dir = 'imagery/'+zone+'/'
if not os.path.exists(img_dir):
	os.mkdir(img_dir)

### Grid dimensions:
# the y-direction
x_start = 0.
y_start = 0.
#x_end = 120000.
#y_end = 40000.

if parameters == 0:
	z = 100

	x_end = 430e3
	y_end = x_end

	nx = 101
	ny = nx

	x_res = x_end/(nx-1)
	y_res = y_end/(ny-1)

elif parameters == 1:
	z = 50

	x_end = 90e3
	x_start = x_end-125e3
	#x_end = 125e3
	y_end = 125e3

	nx = 101
	ny = nx

	x_res = x_end/(nx-1)
	y_res = y_end/(ny-1)

print('x: Number of points:'+str(nx)+' / Resolution: '+str(x_res)+' m')
print('y: Number of points:'+str(ny)+' / Resolution: '+str(y_res)+' m')

# Dimension vectors
x_vec = np.linspace(np.int_(x_start),np.int_(x_end),np.int_(nx))
y_vec = np.linspace(np.int_(y_start),np.int_(y_end),np.int_(ny))
#x_vec = np.arange(x_start, x_end+x_res, x_res)
#y_vec = np.arange(y_start, y_end+y_res, y_res)

# Create the grid
(x_grid, y_grid) = np.meshgrid(x_vec, y_vec)

day2sec = 3600*24
#days_vec = np.array([0, 5, 10, 15, 20, 30, 50, 100])
days_vec = np.arange(12)
#days_vec = [0]
time_vec = day2sec*days_vec

LAT = 38.






###--------------------
### Tang's model run
###--------------------

print('Calculation time set to: '+str(days_vec[len(days_vec)-1])+' days')

# Initialize variables:	
psi     = np.empty([len(x_vec),len(y_vec),len(time_vec)])
psi[:]  = np.nan
ug      = np.empty([len(x_vec),len(y_vec),len(time_vec)])
ug[:]   = np.nan
vg      = np.empty([len(x_vec),len(y_vec),len(time_vec)])
vg[:]   = np.nan
ua      = np.empty([len(x_vec),len(y_vec),len(time_vec)])
ua[:]   = np.nan
va      = np.empty([len(x_vec),len(y_vec),len(time_vec)])
va[:]   = np.nan
w       = np.empty([len(x_vec),len(y_vec),len(time_vec)])
w[:]    = np.nan
div     = np.empty([len(x_vec),len(y_vec),len(time_vec)])
div[:]  = np.nan
vort    = np.empty([len(x_vec),len(y_grid),len(time_vec)])
vort[:] = np.nan

# Loop through the time set
for indt in range(0,len(time_vec)):
	time = time_vec[indt]
	day = days_vec[indt]
	
	# Calculate the Stream function and all the derivated variables for our grid:
	(psi[:,:,indt],ug[:,:,indt],vg[:,:,indt],ua[:,:,indt],va[:,:,indt],w[:,:,indt],div[:,:,indt],vort[:,:,indt]) = tm_equations(x_grid,y_grid,z,time,parameters)

# Velocity fields
ut = ug+ua
vt = vg+va

Ugt = np.sqrt(ug*ug+vg*vg)
Uat = np.sqrt(ua*ua+va*va)
U = np.sqrt(ut*ut+vt*vt)

# Divergence and Vorticity normalization
OMEGA = 7.3e-5
Rt = 6371e3
f = 2*OMEGA*np.sin(np.radians(LAT))
div = div/f
vort = vort/f



# To get the correct Dynamic Height field (in cm)
dh = psi*1e-4/9.8*100
maxvalue = 24
adjustment = maxvalue-maxvalue*y_grid/y_end
dh_adj = np.empty(dh.shape)

w_max = 0
div_max = 0
vort_max = 0
for indt in range(0,len(time_vec)):
	dh_adj[:,:,indt] = dh[:,:,indt] + adjustment
	print(' ')
	print('Day: '+str("{0:02d}".format(days_vec[indt])))
	print('Ug_max = '+str("{:02.3f}".format(np.amax(abs(Ugt[:,:,indt]))))+' || Ua_max = '+str("{:02.3f}".format(np.amax(abs(Uat[:,:,indt]))))+' || w_max(m/day) = '+str(np.amax(abs(w[:,:,indt]))*day2sec))
	print('ua_max = '+str("{:02.3f}".format(np.amax(abs(ua[:,:,indt]))))+' || va_max = '+str(np.amax(abs(va[:,:,indt]))))
	print('div_max/f = '+str(np.amax(div[:,:,indt]))+' || vort_max/f = '+str(np.amax(abs(vort[:,:,indt]))))
	w_max =  max(w_max,np.amax(abs(w[:,:,indt]))*day2sec)
	div_max = max(div_max,np.amax(abs(div[:,:,indt])))
	vort_max = max(vort_max,np.amax(abs(vort[:,:,indt])))
#dh_adj = dh+0.6*y_grid*1e-4/9.8*100-200

# udiff = np.empty(psi.shape)
# vdiff = np.empty(psi.shape)
# for ix in range(x_grid.shape[0]-1):
# 	for iy in range(x_grid.shape[1]-1):
# 		udiff[ix,iy,:] = -(psi[ix+1,iy] - psi[ix,iy])/(y_vec[ix+1] - y_vec[ix])		
# 		vdiff[ix,iy,:] = (psi[ix,iy+1] - psi[ix,iy])/(x_vec[iy+1] - x_vec[iy])












###-------------------
### Ploting Section
###-------------------

dhl = 15
dhn = 101


if plot_dh:
	pl = 'dh'
	print('Ploting '+pl)
	# Create dir
	image_dir = img_dir+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	lev  = np.linspace(-dhl,dhl,dhn)
	colormap = plt.cm.RdYlBu_r
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]	
		pf = dh[:,:,indt]	
		tm_cfplot(x_grid, y_grid, day, plotfield=pf, levels=lev, plotlabel=pl,units='cm',cmap=colormap,image_dir=image_dir)
	#Create GIF
	if plot_gif:
		print('Creating '+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=pl)


if plot_vec_Ugt:
	pl = 'Ugt'
	print('Ploting '+pl)
	# Create dir
	image_dir = img_dir+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	r = 5
	l = 0.5  # length of legend vector
	s = 10   # scale of vector
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]
		vecfield = np.array([ug[:,:,indt],vg[:,:,indt]])
		tm_vecplot(x_grid,y_grid,day,vecfield,veclabel=pl,r=r,length=l,scale=s,units='m/s',image_dir=image_dir)
	#Create GIF
	if plot_gif:
		print('Creating '+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=pl)


if plot_vec_Uat:
	pl = 'Uat'
	print('Ploting '+pl)
	# Create dir
	image_dir = img_dir+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	r = 5
	l = 0.05  # length of legend vector
	s = 1   # scale of vector
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]
		vecfield = np.array([ua[:,:,indt],va[:,:,indt]])
		tm_vecplot(x_grid,y_grid,day,vecfield,veclabel=pl,r=r,length=l,scale=s,units='m/s',image_dir=image_dir)
	#Create GIF
	if plot_gif:
		print('Creating '+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=pl)


if plot_vec_Ut:
	pl = 'Ut'
	print('Ploting '+pl)
	# Create dir
	image_dir = img_dir+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	r = 5
	l = 0.5  # length of legend vector
	s = 10   # scale of vector
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]
		vecfield = np.array([ut[:,:,indt],vt[:,:,indt]])
		tm_vecplot(x_grid,y_grid,day,vecfield,veclabel=pl,r=r,length=l,scale=s,units='m/s',image_dir=image_dir)
	#Create GIF
	if plot_gif:
		print('Creating '+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=pl)


if plot_Ugt_dh:
	vl = 'Ugt'
	pl = 'dh'
	print('Ploting '+vl+'+'+pl)
	# Create dir
	image_dir = img_dir+vl+'+'+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	r = 5
	l = 0.5  # length of legend vector
	s = 5   # scale of vector
	lev  = np.linspace(-dhl,dhl,dhn)
	colormap = plt.cm.RdYlBu_r
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]
		vecfield = np.array([ug[:,:,indt],vg[:,:,indt]])
		pf = dh[:,:,indt]
		tm_vec_cfplot(x_grid,y_grid,day,vecfield,veclabel=vl,vecunits='m/s',r=r,length=l,scale=s,plotfield=pf,plotlabel=pl,levels=lev,plotunits='cm',cmap=colormap,image_dir=image_dir)
	#Create GIF
	if plot_gif:
		print('Creating '+vl+'+'+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=str(vl+'+'+pl))


if plot_Uat_dh:
	vl = 'Uat'
	pl = 'dh'
	print('Ploting '+vl+'+'+pl)
	# Create dir
	image_dir = img_dir+vl+'+'+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	r = 5
	l = 0.01  # length of legend vector
	s = 0.2   # scale of vector
	lev  = np.linspace(-dhl,dhl,dhn)
	colormap = plt.cm.RdYlBu_r
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]
		vecfield = np.array([ua[:,:,indt],va[:,:,indt]])
		pf = dh[:,:,indt]
		tm_vec_cfplot(x_grid,y_grid,day,vecfield,veclabel=vl,vecunits='m/s',r=r,length=l,scale=s,plotfield=pf,plotlabel=pl,levels=lev,plotunits='cm',cmap=colormap,image_dir=image_dir)
	#Create GIF
	if plot_gif:
		print('Creating '+vl+'+'+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=str(vl+'+'+pl))


if plot_Ut_dh:
	vl = 'Ut'
	pl = 'dh'
	print('Ploting '+vl+'+'+pl)
	# Create dir
	image_dir = img_dir+vl+'+'+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	r = 5
	l = 0.5  # length of legend vector
	s = 5   # scale of vector
	lev  = np.linspace(-dhl,dhl,dhn)
	colormap = plt.cm.RdYlBu_r
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]
		vecfield = np.array([ut[:,:,indt],vt[:,:,indt]])
		pf = dh[:,:,indt]
		tm_vec_cfplot(x_grid,y_grid,day,vecfield,veclabel=vl,vecunits='m/s',r=r,length=l,scale=s,plotfield=pf,plotlabel=pl,levels=lev,plotunits='cm',cmap=colormap,image_dir=image_dir)
	#Create GIF
	if plot_gif:
		print('Creating '+vl+'+'+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=str(vl+'+'+pl))


if plot_div_dh:
	pl = 'div+dh'
	print('Ploting '+pl)
	# Create dir
	image_dir = img_dir+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	#levels  = np.linspace(-div_max,div_max,101)
	levels = np.linspace(-3.1,3.1,101)

	cmap = plt.cm.coolwarm
	plotlabel = pl
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]	
		plotfield = div[:,:,indt]*1e3
		plotfield2 = dh[:,:,indt]

		# Plot
		plt.cla
		plt.clf
		plt.close()
		plt.figure()
		# Filled Contour plot
		cs = plt.contourf(x_grid/1000, y_grid/1000, plotfield, levels=levels ,cmap=cmap)
		# Contour lines on top
		cs2= plt.contour(x_grid/1000, y_grid/1000, plotfield2, colors='grey', linewidths=1)
		#plt.clabel(cs2, inline=True, fmt="%0.3f", fontsize=10, colors="black")
		# Colorbar
		#cbar = plt.colorbar(cs, shrink=1.0, extend='both',format='%.0e')
		cbar = plt.colorbar(cs, shrink=1.0, extend='both',ticks=np.int_(np.linspace(-3,3,7)))
		cbar.ax.set_ylabel('div/f (x10$^{-3}$)')

		plt.title('Model '+plotlabel+' d'+str(day))
		plt.ylabel('y (km)')
		plt.xlabel('x (km)')

		image_name = 'model_'+plotlabel+'_d'+str("{0:02d}".format(day))
		filename = image_dir+image_name+'.pdf'
		plt.savefig(filename,bbox_inches='tight')
		filename = image_dir+image_name+'.png'
		plt.savefig(filename,bbox_inches='tight')
		plt.close()

	#Create GIF
	if plot_gif:
		print('Creating '+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=pl)


if plot_Uat_div:
	vl = 'Uat'
	pl = 'div'
	print('Ploting '+vl+'+'+pl)
	# Create dir
	image_dir = img_dir+vl+'+'+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	r = 5
	l = 0.05  # length of legend vector
	s = 1   # scale of vector
	lev  = np.linspace(-div_max,div_max,101)
	colormap = plt.cm.bwr
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]
		vecfield = np.array([ua[:,:,indt],va[:,:,indt]])
		pf = div[:,:,indt]
		tm_vec_cfplot(x_grid,y_grid,day,vecfield,veclabel=vl,vecunits='m/s',r=r,length=l,scale=s,plotfield=pf,plotlabel=pl,levels=lev,plotunits='s-1',cmap=colormap,image_dir=image_dir)
	#Create GIF
	if plot_gif:
		print('Creating '+vl+'+'+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=str(vl+'+'+pl))


if plot_Uat_div_dh:
	vl = 'Uat'
	pl = 'div+dh'
	print('Ploting '+vl+'+'+pl)
	# Create dir
	image_dir = img_dir+vl+'+'+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	r = 5
	length = 0.01  # length of legend vector
	scale = 0.2   # scale of vector
	#levels = np.linspace(-div_max,div_max,101)
	levels = np.linspace(-3.1,3.1,101)
	cmap = plt.cm.coolwarm
	plotlabel = pl
	veclabel = vl
	vecunits = 'm/s'
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]
		vecfield = np.array([ua[:,:,indt],va[:,:,indt]])
		plotfield = div[:,:,indt]*1e3
		plotfield2 = dh[:,:,indt]
		
		# Plot
		plt.cla
		plt.clf
		plt.close()
		plt.figure()

		# Filled Contour plot
		cs = plt.contourf(x_grid/1000, y_grid/1000, plotfield, levels=levels ,cmap=cmap)
		cbar = plt.colorbar(cs, shrink=1.0, extend='both',ticks=np.int_(np.linspace(-3,3,7)))
		cbar.ax.set_ylabel('div/f (x10$^{-3}$)')
		# Contour lines on top
		cs2= plt.contour(x_grid/1000, y_grid/1000, plotfield2, colors='grey', linewidths=1)

		# Vector plot:
		# To obtain a clearer vector map will plot vectors each r data point:
		Q = plt.quiver(x_grid[::r,::r][0,:]/1000,y_grid[::r,::r][:,0]/1000,vecfield[0][::r,::r],vecfield[1][::r,::r],scale=scale)
		qk = plt.quiverkey(Q, 0.65, 0.92, length, str(length)+vecunits, labelpos='E',coordinates='figure')

		plt.title('Model '+veclabel+'+'+plotlabel+' d'+str("{0:02d}".format(day)))
		plt.ylabel('y (km)')
		plt.xlabel('x (km)')

		image_name = 'model_'+veclabel+'+'+plotlabel+'_d'+str("{0:02d}".format(day))+'.pdf'
		filename = image_dir+image_name
		plt.savefig(filename,bbox_inches='tight')
	  	image_name = 'model_'+veclabel+'+'+plotlabel+'_d'+str("{0:02d}".format(day))+'.png'
		filename = image_dir+image_name
		plt.savefig(filename,bbox_inches='tight')

		plt.close()
	#Create GIF
	if plot_gif:
		print('Creating '+vl+'+'+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=str(vl+'+'+pl))


if plot_div_w:
	pl = 'div+w'
	print('Ploting '+pl)
	# Create dir
	image_dir = img_dir+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	levels  = np.linspace(-div_max,div_max,41)
	cmap = plt.cm.bwr
	plotlabel = pl
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]	
		plotfield = div[:,:,indt]
		plotfield2 = w[:,:,indt]*day2sec

		# Plot
		plt.cla
		plt.clf
		plt.close()
		plt.figure()
		# Filled Contour plot
		cs = plt.contourf(x_grid/1000, y_grid/1000, plotfield, levels=levels ,cmap=cmap)
		# Contour lines on top
		cs2= plt.contour(x_grid/1000, y_grid/1000, plotfield2, levels=np.linspace(-w_max,w_max,11),colors='grey', linewidths=1)
		#plt.clabel(cs2, inline=True, fmt="%0.3f", fontsize=10, colors="black")
		# Colorbar
		cbar = plt.colorbar(cs, shrink=1.0, extend='both')
		cbar.ax.set_ylabel('div/f')

		plt.title('Model '+plotlabel+' d'+str(day))
		plt.ylabel('y (km)')
		plt.xlabel('x (km)')

		image_name = 'model_'+plotlabel+'_d'+str("{0:02d}".format(day))
		filename = image_dir+image_name+'.pdf'
		plt.savefig(filename,bbox_inches='tight')
		filename = image_dir+image_name+'.png'
		plt.savefig(filename,bbox_inches='tight')
		plt.close()

	#Create GIF
	if plot_gif:
		print('Creating '+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=pl)


if plot_w_dh:
	pl = 'w+dh'
	print('Ploting '+pl)
	# Create dir
	image_dir = img_dir+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	levels  = np.linspace(-40,40,101)
	cmap = plt.cm.bwr
	plotlabel = pl
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]	
		plotfield2 = dh[:,:,indt]
		plotfield = w[:,:,indt]*day2sec

		# Plot
		plt.cla
		plt.clf
		plt.close()
		plt.figure()
		# Filled Contour plot
		cs = plt.contourf(x_grid/1000, y_grid/1000, plotfield, levels=levels ,cmap=cmap)
		# Contour lines on top
		cs2= plt.contour(x_grid/1000, y_grid/1000, plotfield2, colors='grey', linewidths=1)
		#plt.clabel(cs2, inline=True, colors="grey")
		# Colorbar
		cbar = plt.colorbar(cs, shrink=1.0, extend='both',ticks=np.linspace(-40,40,9))
		cbar.ax.set_ylabel('w (m/day)')

		plt.title('Model '+plotlabel+' d'+str(day))
		plt.ylabel('y (km)')
		plt.xlabel('x (km)')

		image_name = 'model_'+plotlabel+'_d'+str("{0:02d}".format(day))
		filename = image_dir+image_name+'.pdf'
		plt.savefig(filename,bbox_inches='tight')
		filename = image_dir+image_name+'.png'
		plt.savefig(filename,bbox_inches='tight')
		plt.close()

	#Create GIF
	if plot_gif:
		print('Creating '+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=pl)


if plot_vort_dh:
	pl = 'vort+dh'
	print('Ploting '+pl)
	# Create dir
	image_dir = img_dir+pl+'/'
	if not os.path.exists(image_dir):
		os.mkdir(image_dir)
	#levels  = np.linspace(-0.16,0.16,101)
	levels = np.linspace(-1.6,1.6,101)
	cmap = plt.cm.coolwarm
	plotlabel = pl
	for indt in range(0,len(time_vec)):
		time = time_vec[indt]
		day = days_vec[indt]	
		plotfield = vort[:,:,indt]*1e1
		plotfield2 = dh[:,:,indt]

		# Plot
		plt.cla
		plt.clf
		plt.close()
		plt.figure()
		# Filled Contour plot
		cs = plt.contourf(x_grid/1000, y_grid/1000, plotfield, levels=levels ,cmap=cmap)
		# Contour lines on top
		cs2= plt.contour(x_grid/1000, y_grid/1000, plotfield2, colors='grey', linewidths=1)
		#plt.clabel(cs2, inline=True, fmt="%0.3f", fontsize=10, colors="black")
		# Colorbar
		cbar = plt.colorbar(cs, shrink=1.0, extend='both',ticks=np.linspace(-1.6,1.6,5))
		cbar.ax.set_ylabel('vort/f (x10$^{-1}$)')

		plt.title('Model '+plotlabel+' d'+str(day))
		plt.ylabel('y (km)')
		plt.xlabel('x (km)')

		image_name = 'model_'+plotlabel+'_d'+str("{0:02d}".format(day))
		filename = image_dir+image_name+'.pdf'
		plt.savefig(filename,bbox_inches='tight')
		filename = image_dir+image_name+'.png'
		plt.savefig(filename,bbox_inches='tight')
		plt.close()

	#Create GIF
	if plot_gif:
		print('Creating '+pl+' GIF')
		tm_gif(days_vec,img_dir,plotlabel=pl)







