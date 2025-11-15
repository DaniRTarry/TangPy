"""
Plotting functions for tangmodel.py

Author: Tarry
Type: FUNCTION
Description: Contains all plotting functions required in tangmodel.py
"""

import matplotlib.pyplot as plt
import numpy as np
import os


# Constants
KM_CONVERSION = 1000  # Conversion from meters to kilometers


def _save_figure(image_dir, image_name, formats=('pdf', 'png')):
    """Helper function to save figures in multiple formats."""
    for fmt in formats:
        filename = os.path.join(image_dir, f"{image_name}.{fmt}")
        plt.savefig(filename, bbox_inches='tight')


def _setup_colorbar(cs, label, units, levels):
    """Helper function to create and configure colorbar."""
    ticks = np.linspace(levels[0], abs(levels[0]), 11)
    cbar = plt.colorbar(cs, shrink=1.0, extend='both', ticks=ticks)
    cbar.ax.set_ylabel(f'{label} ({units})')
    return cbar


def tm_cfplot(x_grid, y_grid, day, plotfield, levels, plotlabel, units, cmap, image_dir):
    """
    Create filled contour plot.
    
    Parameters:
        x_grid: X coordinate grid
        y_grid: Y coordinate grid
        day: Day number for plot title
        plotfield: Data field to plot
        levels: Contour levels
        plotlabel: Label for the plot
        units: Units for the data
        cmap: Colormap
        image_dir: Directory to save images
    """
    plt.close('all')
    plt.figure()
    
    # Filled contour plot
    cs = plt.contourf(x_grid / KM_CONVERSION, y_grid / KM_CONVERSION, 
                      plotfield, levels=levels, cmap=cmap)
    
    # Colorbar
    _setup_colorbar(cs, plotlabel, units, levels)
    
    # Labels and title
    plt.title(f'Model {plotlabel} d{day}')
    plt.ylabel('y (km)')
    plt.xlabel('x (km)')
    
    # Save figure
    image_name = f'model_{plotlabel}_d{day:02d}'
    _save_figure(image_dir, image_name)
    plt.close()


def tm_vec_cfplot(x_grid, y_grid, day, vecfield, veclabel, vecunits, r, 
                  length, scale, plotfield, plotlabel, levels, plotunits, 
                  cmap, image_dir):
    """
    Create combined vector and filled contour plot.
    
    Parameters:
        x_grid: X coordinate grid
        y_grid: Y coordinate grid
        day: Day number for plot title
        vecfield: Vector field data [u, v]
        veclabel: Label for vector field
        vecunits: Units for vector field
        r: Stride for vector plotting (plot every r-th point)
        length: Reference vector length
        scale: Scale for quiver plot
        plotfield: Background contour field
        plotlabel: Label for contour field
        levels: Contour levels
        plotunits: Units for contour field
        cmap: Colormap
        image_dir: Directory to save images
    """
    plt.close('all')
    plt.figure()
    
    # Filled contour plot
    cs = plt.contourf(x_grid / KM_CONVERSION, y_grid / KM_CONVERSION, 
                      plotfield, levels=levels, cmap=cmap)
    _setup_colorbar(cs, plotlabel, plotunits, levels)
    
    # Vector plot (subsample for clarity)
    Q = plt.quiver(x_grid[::r, ::r][0, :] / KM_CONVERSION,
                   y_grid[::r, ::r][:, 0] / KM_CONVERSION,
                   vecfield[0][::r, ::r],
                   vecfield[1][::r, ::r],
                   scale=scale, headwidth=4, headlength=4)
    
    plt.quiverkey(Q, 0.65, 0.92, length, f'{length}{vecunits}', 
                  labelpos='E', coordinates='figure')
    
    # Labels and title
    plt.title(f'Model {veclabel}+{plotlabel} d{day:02d}')
    plt.ylabel('y (km)')
    plt.xlabel('x (km)')
    
    # Save figure
    image_name = f'model_{veclabel}+{plotlabel}_d{day:02d}'
    _save_figure(image_dir, image_name)
    plt.close()


def tm_vecplot(x_grid, y_grid, day, vecfield, veclabel, r, length, 
               scale, units, image_dir):
    """
    Create vector plot.
    
    Parameters:
        x_grid: X coordinate grid
        y_grid: Y coordinate grid
        day: Day number for plot title
        vecfield: Vector field data [u, v]
        veclabel: Label for vector field
        r: Stride for vector plotting (plot every r-th point)
        length: Reference vector length
        scale: Scale for quiver plot
        units: Units for vector field
        image_dir: Directory to save images
    """
    plt.close('all')
    plt.figure()
    
    # Vector plot (subsample for clarity)
    Q = plt.quiver(x_grid[::r, ::r] / KM_CONVERSION,
                   y_grid[::r, ::r] / KM_CONVERSION,
                   vecfield[0][::r, ::r],
                   vecfield[1][::r, ::r],
                   scale=scale)
    
    plt.quiverkey(Q, 0.9, 0.95, length, f'{length}{units}', 
                  labelpos='E', coordinates='figure')
    
    # Labels and title
    plt.title(f'Model {veclabel} d{day:02d}')
    plt.ylabel('y (km)')
    plt.xlabel('x (km)')
    
    # Save figure
    image_name = f'model_{veclabel}_d{day:02d}'
    _save_figure(image_dir, image_name)
    plt.close()
