"""
Main program for Tang's model field generation

Author: Tarry
Type: PROGRAM
Description: Generates and visualizes fields from Tang's model
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from typing import Dict, List, Tuple

from tang_model_equations import tm_equations
from tang_model_plot import tm_cfplot, tm_vec_cfplot, tm_vecplot
from gifcreator import tm_gif


# Physical constants
OMEGA_EARTH = 7.3e-5
EARTH_RADIUS = 6371e3
LATITUDE = 38.0
DAY_TO_SECONDS = 3600 * 24


# Configuration
PLOT_CONFIG = {
    'psi': False,
    'dh': True,
    'vec_Ugt': False,
    'vec_Uat': False,
    'vec_Ut': False,
    'Ugt_dh': True,
    'Uat_dh': True,
    'Ut_dh': False,
    'div_dh': True,
    'Uat_div': False,
    'Uat_div_dh': True,
    'div_w': False,
    'w_dh': True,
    'vort_dh': True,
    'generate_gifs': True,
}

# Model parameters
PARAMETER_CONFIGS = {
    0: {
        'name': 'Tang_paper',
        'z': 100,
        'x_start': 0.0,
        'x_end': 430e3,
        'y_start': 0.0,
        'y_end': 430e3,
        'nx': 101,
        'ny': 101,
    },
    1: {
        'name': 'Algerian_Basin',
        'z': 50,
        'x_start': -35e3,
        'x_end': 90e3,
        'y_start': 0.0,
        'y_end': 125e3,
        'nx': 101,
        'ny': 101,
    },
}

# Plot settings for different field types
PLOT_SETTINGS = {
    'dh': {
        'levels': (-15, 15, 101),
        'colormap': 'RdYlBu_r',
        'units': 'cm',
    },
    'div': {
        'levels': (-3.1, 3.1, 101),
        'colormap': 'coolwarm',
        'units': 's-1',
        'label': 'div/f (x10$^{-3}$)',
        'scale_factor': 1e3,
    },
    'w': {
        'levels': (-40, 40, 101),
        'colormap': 'bwr',
        'units': 'm/day',
        'scale_factor': DAY_TO_SECONDS,
    },
    'vort': {
        'levels': (-1.6, 1.6, 101),
        'colormap': 'coolwarm',
        'units': 'f-1',
        'label': 'vort/f (x10$^{-1}$)',
        'scale_factor': 1e1,
    },
}

# Vector plot settings
VECTOR_SETTINGS = {
    'Ugt': {'r': 5, 'length': 0.5, 'scale': 10, 'units': 'm/s'},
    'Uat': {'r': 5, 'length': 0.05, 'scale': 1, 'units': 'm/s'},
    'Ut': {'r': 5, 'length': 0.5, 'scale': 10, 'units': 'm/s'},
}


def setup_grid(config: Dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Create spatial grid for model domain.
    
    Returns:
        Tuple of (x_vec, y_vec, x_grid, y_grid)
    """
    x_vec = np.linspace(config['x_start'], config['x_end'], config['nx'])
    y_vec = np.linspace(config['y_start'], config['y_end'], config['ny'])
    x_grid, y_grid = np.meshgrid(x_vec, y_vec)
    
    x_res = (config['x_end'] - config['x_start']) / (config['nx'] - 1)
    y_res = (config['y_end'] - config['y_start']) / (config['ny'] - 1)
    
    print(f"x: Number of points: {config['nx']} / Resolution: {x_res:.2f} m")
    print(f"y: Number of points: {config['ny']} / Resolution: {y_res:.2f} m")
    
    return x_vec, y_vec, x_grid, y_grid


def initialize_fields(nx: int, ny: int, nt: int) -> Dict[str, np.ndarray]:
    """Initialize all field arrays with NaN."""
    shape = (nx, ny, nt)
    fields = {}
    for field in ['psi', 'ug', 'vg', 'ua', 'va', 'w', 'div', 'vort']:
        fields[field] = np.full(shape, np.nan)
    return fields


def calculate_fields(x_grid: np.ndarray, y_grid: np.ndarray, z: float,
                    time_vec: np.ndarray, parameter_id: int) -> Dict[str, np.ndarray]:
    """
    Calculate all model fields for specified times.
    
    Returns:
        Dictionary containing all calculated fields
    """
    nx, ny = x_grid.shape
    nt = len(time_vec)
    fields = initialize_fields(nx, ny, nt)
    
    # Calculate fields at each time step
    for indt, time in enumerate(time_vec):
        results = tm_equations(x_grid, y_grid, z, time, parameter_id)
        fields['psi'][:, :, indt] = results[0]
        fields['ug'][:, :, indt] = results[1]
        fields['vg'][:, :, indt] = results[2]
        fields['ua'][:, :, indt] = results[3]
        fields['va'][:, :, indt] = results[4]
        fields['w'][:, :, indt] = results[5]
        fields['div'][:, :, indt] = results[6]
        fields['vort'][:, :, indt] = results[7]
    
    return fields


def calculate_derived_fields(fields: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Calculate derived velocity and scalar fields."""
    derived = {}
    
    # Total velocities
    derived['ut'] = fields['ug'] + fields['ua']
    derived['vt'] = fields['vg'] + fields['va']
    
    # Velocity magnitudes
    derived['Ugt'] = np.sqrt(fields['ug']**2 + fields['vg']**2)
    derived['Uat'] = np.sqrt(fields['ua']**2 + fields['va']**2)
    derived['U'] = np.sqrt(derived['ut']**2 + derived['vt']**2)
    
    # Normalize divergence and vorticity by Coriolis parameter
    f = 2 * OMEGA_EARTH * np.sin(np.radians(LATITUDE))
    derived['div_f'] = fields['div'] / f
    derived['vort_f'] = fields['vort'] / f
    
    # Dynamic height (in cm)
    derived['dh'] = fields['psi'] * 1e-4 / 9.8 * 100
    
    return derived


def print_field_statistics(derived: Dict[str, np.ndarray], days_vec: np.ndarray):
    """Print maximum values for key fields at each time step."""
    w_max = div_max = vort_max = 0
    
    for indt, day in enumerate(days_vec):
        print(f'\nDay: {day:02d}')
        print(f"Ug_max = {np.amax(np.abs(derived['Ugt'][:,:,indt])):.3f} || "
              f"Ua_max = {np.amax(np.abs(derived['Uat'][:,:,indt])):.3f} || "
              f"w_max(m/day) = {np.amax(np.abs(derived['w'][:,:,indt])) * DAY_TO_SECONDS:.3f}")
        print(f"ua_max = {np.amax(np.abs(derived['ua'][:,:,indt])):.3f} || "
              f"va_max = {np.amax(np.abs(derived['va'][:,:,indt])):.3f}")
        print(f"div_max/f = {np.amax(derived['div_f'][:,:,indt]):.6f} || "
              f"vort_max/f = {np.amax(np.abs(derived['vort_f'][:,:,indt])):.6f}")
        
        w_max = max(w_max, np.amax(np.abs(derived['w'][:,:,indt])) * DAY_TO_SECONDS)
        div_max = max(div_max, np.amax(np.abs(derived['div_f'][:,:,indt])))
        vort_max = max(vort_max, np.amax(np.abs(derived['vort_f'][:,:,indt])))
    
    return w_max, div_max, vort_max


def ensure_directory(path: str):
    """Create directory if it doesn't exist."""
    Path(path).mkdir(parents=True, exist_ok=True)


def plot_contour_field(x_grid: np.ndarray, y_grid: np.ndarray, days_vec: np.ndarray,
                       field: np.ndarray, field_name: str, settings: Dict,
                       img_dir: str, create_gif: bool = True):
    """Generic function to plot contour fields."""
    print(f'Plotting {field_name}')
    image_dir = os.path.join(img_dir, field_name, '')
    ensure_directory(image_dir)
    
    levels = np.linspace(*settings['levels'])
    colormap = getattr(plt.cm, settings['colormap'])
    scale = settings.get('scale_factor', 1.0)
    
    for indt, day in enumerate(days_vec):
        plot_data = field[:, :, indt] * scale
        tm_cfplot(x_grid, y_grid, day, plotfield=plot_data, levels=levels,
                 plotlabel=field_name, units=settings['units'],
                 cmap=colormap, image_dir=image_dir)
    
    if create_gif:
        print(f'Creating {field_name} GIF')
        tm_gif(days_vec, img_dir, plotlabel=field_name)


def plot_vector_field(x_grid: np.ndarray, y_grid: np.ndarray, days_vec: np.ndarray,
                     u_field: np.ndarray, v_field: np.ndarray, field_name: str,
                     settings: Dict, img_dir: str, create_gif: bool = True):
    """Generic function to plot vector fields."""
    print(f'Plotting {field_name}')
    image_dir = os.path.join(img_dir, field_name, '')
    ensure_directory(image_dir)
    
    for indt, day in enumerate(days_vec):
        vecfield = np.array([u_field[:, :, indt], v_field[:, :, indt]])
        tm_vecplot(x_grid, y_grid, day, vecfield, veclabel=field_name,
                  r=settings['r'], length=settings['length'],
                  scale=settings['scale'], units=settings['units'],
                  image_dir=image_dir)
    
    if create_gif:
        print(f'Creating {field_name} GIF')
        tm_gif(days_vec, img_dir, plotlabel=field_name)


def plot_vector_contour_combined(x_grid: np.ndarray, y_grid: np.ndarray,
                                 days_vec: np.ndarray, u_field: np.ndarray,
                                 v_field: np.ndarray, contour_field: np.ndarray,
                                 vec_name: str, contour_name: str,
                                 vec_settings: Dict, contour_settings: Dict,
                                 img_dir: str, create_gif: bool = True):
    """Plot vector field overlaid on contour field."""
    combined_name = f'{vec_name}+{contour_name}'
    print(f'Plotting {combined_name}')
    image_dir = os.path.join(img_dir, combined_name, '')
    ensure_directory(image_dir)
    
    levels = np.linspace(*contour_settings['levels'])
    colormap = getattr(plt.cm, contour_settings['colormap'])
    scale = contour_settings.get('scale_factor', 1.0)
    
    for indt, day in enumerate(days_vec):
        vecfield = np.array([u_field[:, :, indt], v_field[:, :, indt]])
        plot_data = contour_field[:, :, indt] * scale
        
        tm_vec_cfplot(x_grid, y_grid, day, vecfield, veclabel=vec_name,
                     vecunits=vec_settings['units'], r=vec_settings['r'],
                     length=vec_settings['length'], scale=vec_settings['scale'],
                     plotfield=plot_data, plotlabel=contour_name,
                     levels=levels, plotunits=contour_settings['units'],
                     cmap=colormap, image_dir=image_dir)
    
    if create_gif:
        print(f'Creating {combined_name} GIF')
        tm_gif(days_vec, img_dir, plotlabel=combined_name)


def plot_special_combined(x_grid: np.ndarray, y_grid: np.ndarray, days_vec: np.ndarray,
                         field1: np.ndarray, field2: np.ndarray, plot_name: str,
                         settings: Dict, img_dir: str, create_gif: bool = True):
    """Plot special combined fields with contour lines (div+dh, w+dh, vort+dh)."""
    print(f'Plotting {plot_name}')
    image_dir = os.path.join(img_dir, plot_name, '')
    ensure_directory(image_dir)
    
    levels = np.linspace(*settings['levels'])
    colormap = getattr(plt.cm, settings['colormap'])
    scale = settings.get('scale_factor', 1.0)
    
    for indt, day in enumerate(days_vec):
        plt.close('all')
        plt.figure()
        
        # Filled contour plot
        plot_data = field1[:, :, indt] * scale
        cs = plt.contourf(x_grid / 1000, y_grid / 1000, plot_data,
                         levels=levels, cmap=colormap)
        
        # Contour lines
        plt.contour(x_grid / 1000, y_grid / 1000, field2[:, :, indt],
                   colors='grey', linewidths=1)
        
        # Colorbar
        if 'label' in settings:
            cbar = plt.colorbar(cs, shrink=1.0, extend='both',
                              ticks=np.int_(np.linspace(levels[0], levels[-1], 7)))
            cbar.ax.set_ylabel(settings['label'])
        else:
            cbar = plt.colorbar(cs, shrink=1.0, extend='both')
            cbar.ax.set_ylabel(f"{plot_name} ({settings['units']})")
        
        plt.title(f'Model {plot_name} d{day:02d}')
        plt.ylabel('y (km)')
        plt.xlabel('x (km)')
        
        image_name = f'model_{plot_name}_d{day:02d}'
        plt.savefig(os.path.join(image_dir, f'{image_name}.pdf'), bbox_inches='tight')
        plt.savefig(os.path.join(image_dir, f'{image_name}.png'), bbox_inches='tight')
        plt.close()
    
    if create_gif:
        print(f'Creating {plot_name} GIF')
        tm_gif(days_vec, img_dir, plotlabel=plot_name)


def main():
    """Main execution function."""
    print('Running file: tang_model.py')
    
    # Select parameter set
    parameter_id = 1
    config = PARAMETER_CONFIGS[parameter_id]
    zone = config['name']
    print(f'Model Parameters set for: {zone}')
    
    # Setup directories
    img_dir = os.path.join('imagery', zone, '')
    ensure_directory(img_dir)
    
    # Setup grid
    x_vec, y_vec, x_grid, y_grid = setup_grid(config)
    
    # Setup time vector
    days_vec = np.arange(12)
    time_vec = DAY_TO_SECONDS * days_vec
    print(f'Calculation time set to: {days_vec[-1]} days')
    
    # Calculate fields
    fields = calculate_fields(x_grid, y_grid, config['z'], time_vec, parameter_id)
    derived = calculate_derived_fields(fields)
    
    # Store fields in derived for convenience
    derived.update(fields)
    
    # Print statistics
    w_max, div_max, vort_max = print_field_statistics(derived, days_vec)
    
    # Generate plots based on configuration
    create_gif = PLOT_CONFIG['generate_gifs']
    
    # Simple contour plots
    if PLOT_CONFIG['dh']:
        plot_contour_field(x_grid, y_grid, days_vec, derived['dh'],
                          'dh', PLOT_SETTINGS['dh'], img_dir, create_gif)
    
    # Vector plots
    if PLOT_CONFIG['vec_Ugt']:
        plot_vector_field(x_grid, y_grid, days_vec, derived['ug'], derived['vg'],
                         'Ugt', VECTOR_SETTINGS['Ugt'], img_dir, create_gif)
    
    if PLOT_CONFIG['vec_Uat']:
        plot_vector_field(x_grid, y_grid, days_vec, derived['ua'], derived['va'],
                         'Uat', VECTOR_SETTINGS['Uat'], img_dir, create_gif)
    
    if PLOT_CONFIG['vec_Ut']:
        plot_vector_field(x_grid, y_grid, days_vec, derived['ut'], derived['vt'],
                         'Ut', VECTOR_SETTINGS['Ut'], img_dir, create_gif)
    
    # Combined vector + contour plots
    if PLOT_CONFIG['Ugt_dh']:
        plot_vector_contour_combined(x_grid, y_grid, days_vec,
                                    derived['ug'], derived['vg'], derived['dh'],
                                    'Ugt', 'dh', VECTOR_SETTINGS['Ugt'],
                                    PLOT_SETTINGS['dh'], img_dir, create_gif)
    
    if PLOT_CONFIG['Uat_dh']:
        vec_settings = VECTOR_SETTINGS['Uat'].copy()
        vec_settings.update({'length': 0.01, 'scale': 0.2})
        plot_vector_contour_combined(x_grid, y_grid, days_vec,
                                    derived['ua'], derived['va'], derived['dh'],
                                    'Uat', 'dh', vec_settings,
                                    PLOT_SETTINGS['dh'], img_dir, create_gif)
    
    if PLOT_CONFIG['Ut_dh']:
        vec_settings = VECTOR_SETTINGS['Ut'].copy()
        vec_settings.update({'scale': 5})
        plot_vector_contour_combined(x_grid, y_grid, days_vec,
                                    derived['ut'], derived['vt'], derived['dh'],
                                    'Ut', 'dh', vec_settings,
                                    PLOT_SETTINGS['dh'], img_dir, create_gif)
    
    # Special combined plots with contour lines
    if PLOT_CONFIG['div_dh']:
        plot_special_combined(x_grid, y_grid, days_vec,
                            derived['div_f'], derived['dh'],
                            'div+dh', PLOT_SETTINGS['div'], img_dir, create_gif)
    
    if PLOT_CONFIG['w_dh']:
        plot_special_combined(x_grid, y_grid, days_vec,
                            derived['w'], derived['dh'],
                            'w+dh', PLOT_SETTINGS['w'], img_dir, create_gif)
    
    if PLOT_CONFIG['vort_dh']:
        plot_special_combined(x_grid, y_grid, days_vec,
                            derived['vort_f'], derived['dh'],
                            'vort+dh', PLOT_SETTINGS['vort'], img_dir, create_gif)
    
    print('\nProcessing complete!')


if __name__ == '__main__':
    main()
