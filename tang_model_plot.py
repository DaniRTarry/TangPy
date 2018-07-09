# -*- coding: utf-8 -*-
"""
@author: Tarry

type: FUNCTION

description: function which contains all the plottings required in tangmodel.py

input: x grid, y grid, day, plotfield, plotlabel
"""

# Load libraries
import matplotlib.pyplot as plt
import numpy as np
import os

# Figures out the absolute path in case the working directory
# my_path = os.path.abspath(__file__)

# A = os.path.join(os.path.dirname(__file__), '..')
# A is the parent directory of the directory where program resides.

# B = os.path.dirname(os.path.realpath(__file__))
# B is the canonicalised (?) directory where the program resides.

# C = os.path.abspath(os.path.dirname(__file__))
# C is the absolute path of the directory where the program resides.


def tm_cfplot(x_grid,y_grid,day,plotfield,levels,plotlabel,units,cmap,image_dir):

	plt.cla
	plt.clf
	plt.close()
	plt.figure()

	# Filled Contour plot
	cs = plt.contourf(x_grid/1000, y_grid/1000, plotfield, levels=levels ,cmap=cmap)
	# Contour lines on top
	#cs2= plt.contour(x_grid/1000, y_grid/1000, plotfield, colors='k', linewidths=1)
	#plt.clabel(cs2, inline=True, fmt="%0.3f", fontsize=10, colors="black")
	# Colorbar
	cbar = plt.colorbar(cs, shrink=1.0, extend='both',ticks=np.linspace(levels[0],abs(levels[0]),11))
	cbar.ax.set_ylabel(plotlabel+' ('+units+')')

	plt.title('Model '+plotlabel+' d'+str(day))
	plt.ylabel('y (km)')
	plt.xlabel('x (km)')

	image_name = 'model_'+plotlabel+'_d'+str("{0:02d}".format(day))
	filename = image_dir+image_name+'.pdf'
	plt.savefig(filename,bbox_inches='tight')
	filename = image_dir+image_name+'.png'
	plt.savefig(filename,bbox_inches='tight')

	plt.close()



def tm_vec_cfplot(x_grid,y_grid,day,vecfield,veclabel,vecunits,r,length,scale,plotfield,plotlabel,levels,plotunits,cmap,image_dir):

	plt.cla
	plt.clf
	plt.close()
	plt.figure()

	# Filled Contour plot
	cs = plt.contourf(x_grid/1000, y_grid/1000, plotfield, levels=levels ,cmap=cmap)
	cbar = plt.colorbar(cs, shrink=1.0, extend='both',ticks=np.linspace(levels[0],abs(levels[0]),11))
	cbar.ax.set_ylabel(plotlabel+' ('+plotunits+')')

	# Vector plot:
	# To obtain a clearer vector map will plot vectors each r data point:
	Q = plt.quiver(x_grid[::r,::r][0,:]/1000,y_grid[::r,::r][:,0]/1000,vecfield[0][::r,::r],vecfield[1][::r,::r],scale=scale,headwidth=4,headlength=4)
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



def tm_vecplot(x_grid,y_grid,day,vecfield,veclabel,r,length,scale,units,image_dir):
	plt.cla
	plt.clf
	plt.close()
	plt.figure()

	# Vector plot:
	# To obtain a clearer vector map will plot vectors each r data point:
	Q = plt.quiver(x_grid[::r,::r]/1000,y_grid[::r,::r]/1000,vecfield[0][::r,::r],vecfield[1][::r,::r],scale=scale)
	qk = plt.quiverkey(Q, 0.9, 0.95, length, str(length)+units, labelpos='E',coordinates='figure')

	plt.title('Model '+veclabel+' d'+str("{0:02d}".format(day)))
	plt.ylabel('y (km)')
	plt.xlabel('x (km)')

	image_name = 'model_'+veclabel+'_d'+str("{0:02d}".format(day))
	filename = image_dir+image_name+'.pdf'
	plt.savefig(filename,bbox_inches='tight')
	filename = image_dir+image_name+'.png'
	plt.savefig(filename,bbox_inches='tight')

	plt.close()





