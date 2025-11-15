"""
Tang Model Physical Parametrization

Author: Tarry
Type: FUNCTION
Description: Contains all physical parametrization of the Tang Model
Input: x1, x2, x3, time, parameters
Output: psi, ug, vg, ua, va, w, divergence, vorticity
"""

import numpy as np
from typing import Tuple, Dict


# Physical constants
OMEGA_EARTH = 7.3e-5  # Earth's angular velocity (rad/s)
EARTH_RADIUS = 6371e3  # Earth's radius (m)
LATITUDE = 38.0  # Latitude (degrees)

# Amplitude parameter
AMPLITUDE = 5000.0


def get_parameter_set(parameter_id: int) -> Dict[str, float]:
    """
    Get physical parameters for the Tang Model.
    
    Parameters:
        parameter_id: 0 or 1, selects parameter configuration
        
    Returns:
        Dictionary containing model parameters
    """
    parameter_sets = {
        0: {
            'wavelength_x': 427e3,      # Wavelength in x1 direction (m)
            'wavelength_y': 427e3,      # Distance between nodal surfaces (m)
            'n1': 6e-3,                 # Brunt-Vaisala frequency layer 1 (rad/s)
            'n2': 2e-3,                 # Brunt-Vaisala frequency layer 2 (rad/s)
            'u0': 0.6,                  # Basic zonal current (m/s)
            'h1': 1e3,                  # Depth of motion layer (m)
            'h2': 1e3,                  # Depth of quiescent layer (m)
        },
        1: {
            'wavelength_x': 250e3,
            'wavelength_y': 125e3,
            'n1': 5.8e-3,
            'n2': 2.3e-6,
            'u0': 0.3,
            'h1': 500.0,
            'h2': 2000.0,
        }
    }
    
    if parameter_id not in parameter_sets:
        raise ValueError(f"Invalid parameter_id: {parameter_id}. Must be 0 or 1.")
    
    return parameter_sets[parameter_id]


def calculate_coriolis_parameters(latitude: float) -> Tuple[float, float]:
    """
    Calculate Coriolis parameter and beta plane approximation.
    
    Parameters:
        latitude: Latitude in degrees
        
    Returns:
        Tuple of (f, beta) where:
            f: Coriolis parameter (rad/s)
            beta: Beta plane parameter (1/(m*s))
    """
    lat_rad = np.radians(latitude)
    f = 2 * OMEGA_EARTH * np.sin(lat_rad)
    beta = 2 * OMEGA_EARTH * np.cos(lat_rad) / EARTH_RADIUS
    return f, beta


def calculate_wave_parameters(params: Dict[str, float], f: float) -> Dict[str, float]:
    """
    Calculate wave numbers and related parameters.
    
    Parameters:
        params: Parameter dictionary from get_parameter_set()
        f: Coriolis parameter
        
    Returns:
        Dictionary of wave parameters
    """
    # Wave numbers
    k = 2.0 * np.pi / params['wavelength_x']
    l = np.pi / params['wavelength_y']
    mu = np.sqrt(l**2 + k**2)
    
    # Stability ratio (wave cutoff at R=3)
    r = params['n1'] / params['n2']
    
    # Vertical wave numbers
    k1 = mu * (params['h1'] * params['n1']) / f
    k2 = mu * (params['h2'] * params['n2']) / f
    
    # Coupling coefficient
    kappa = np.tanh(k2) / np.tanh(k1)
    
    # Auxiliary coefficients
    tanh_k1_half = np.tanh(k1 / 2)
    denom = k1 - np.tanh(k1)
    a = 2 * ((k1 / 2) - tanh_k1_half) / denom
    b = 2 * ((k1 / 2) - 1 / tanh_k1_half) / denom
    
    return {
        'k': k,
        'l': l,
        'mu': mu,
        'r': r,
        'k1': k1,
        'k2': k2,
        'kappa': kappa,
        'a': a,
        'b': b
    }


def calculate_phase_speed(params: Dict[str, float], wave_params: Dict[str, float]) -> Tuple[float, float]:
    """
    Calculate complex phase speed c = cr + i*ci.
    
    Parameters:
        params: Physical parameters
        wave_params: Wave parameters from calculate_wave_parameters()
        
    Returns:
        Tuple of (cr, ci) - real and imaginary parts of phase speed
    """
    u0 = params['u0']
    r = wave_params['r']
    kappa = wave_params['kappa']
    k1 = wave_params['k1']
    a = wave_params['a']
    b = wave_params['b']
    
    tanh_k1 = np.tanh(k1)
    r_kappa = r * kappa
    
    # Real part of phase speed
    cr = (u0 / 2) * (1 - r_kappa * tanh_k1 / (k1 * (1 + r_kappa)))
    
    # Imaginary part of phase speed
    numerator = u0 * (k1 - tanh_k1) * np.sqrt(-(r_kappa + a) * (r_kappa + b))
    denominator = 2 * k1 * (1 + r_kappa)
    ci = numerator / denominator
    
    return cr, ci


def calculate_streamfunction_amplitude(x3: np.ndarray, params: Dict[str, float], 
                                       wave_params: Dict[str, float], 
                                       cr: float, ci: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate real and imaginary parts of streamfunction amplitude.
    
    Returns:
        Tuple of (psi_amplitude, theta) where theta is the phase angle
    """
    u0 = params['u0']
    r = wave_params['r']
    kappa = wave_params['kappa']
    k1 = wave_params['k1']
    h1 = params['h1']
    
    xi = k1 * x3 / h1
    tanh_k1 = np.tanh(k1)
    r_kappa_plus_1 = 1 + r * kappa
    k1_factor = k1 / tanh_k1 - 1
    
    # Real amplitude
    cosh_xi = np.cosh(xi)
    sinh_xi = np.sinh(xi)
    psi_real = AMPLITUDE * (cosh_xi + k1 * ((cr / u0) * r_kappa_plus_1 - 1) * sinh_xi / k1_factor)
    
    # Imaginary amplitude
    psi_imag = AMPLITUDE * k1 * (ci / u0) * r_kappa_plus_1 * sinh_xi / k1_factor
    
    # Amplitude and phase
    psi_amplitude = np.sqrt(psi_real**2 + psi_imag**2)
    theta = np.arctan2(psi_imag, psi_real)
    
    return psi_amplitude, theta


def calculate_vertical_velocity_amplitude(x3: np.ndarray, params: Dict[str, float],
                                          wave_params: Dict[str, float],
                                          cr: float, ci: float, f: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate vertical velocity amplitude and phase.
    
    Returns:
        Tuple of (w_amplitude, sigma) where sigma is the phase angle
    """
    u0 = params['u0']
    n1 = params['n1']
    h1 = params['h1']
    r = wave_params['r']
    kappa = wave_params['kappa']
    k1 = wave_params['k1']
    k = wave_params['k']
    
    xi = k1 * x3 / h1
    r_kappa_plus_1 = 1 + r * kappa
    k1_factor = k1 / np.tanh(k1) - 1
    
    # Auxiliary coefficients
    b_real = k1 * ((cr / u0) * r_kappa_plus_1 - 1) / k1_factor
    b_imag = k1 * (ci / u0) * r_kappa_plus_1 / k1_factor
    
    cosh_xi = np.cosh(xi)
    sinh_xi = np.sinh(xi)
    
    # Real part
    w_real_1 = (k1 * ci / u0 + b_imag) * sinh_xi
    w_real_2 = (xi * b_imag - k1 / u0 * (cr * b_imag + ci * b_real)) * cosh_xi
    w_real = -f * u0 * k * AMPLITUDE / (n1**2 * h1) * (w_real_1 - w_real_2)
    
    # Imaginary part
    w_imag_1 = (xi - k1 * cr / u0 - b_real) * sinh_xi
    w_imag_2 = (xi * b_real - 1 - k1 / u0 * (cr * b_real - ci * b_imag)) * cosh_xi
    w_imag = -f * u0 * k * AMPLITUDE / (n1**2 * h1) * (w_imag_1 + w_imag_2)
    
    # Amplitude and phase
    w_amplitude = np.sqrt(w_real**2 + w_imag**2)
    sigma = np.arctan2(-w_imag, w_real)
    
    return w_amplitude, sigma


def tm_equations(x1: np.ndarray, x2: np.ndarray, x3: np.ndarray, 
                 t: float, parameter_id: int) -> Tuple:
    """
    Calculate Tang Model equations for atmospheric wave dynamics.
    
    Parameters:
        x1: Zonal coordinate array (m)
        x2: Meridional coordinate array (m)
        x3: Vertical coordinate array (m)
        t: Time (s)
        parameter_id: Parameter set selector (0 or 1)
        
    Returns:
        Tuple of (psi, ug, vg, ua, va, w, divergence, vorticity)
            psi: Streamfunction
            ug, vg: Geostrophic velocities
            ua, va: Ageostrophic velocities
            w: Vertical velocity
            divergence: Horizontal divergence
            vorticity: Relative vorticity
    """
    # Get parameters
    params = get_parameter_set(parameter_id)
    f, beta = calculate_coriolis_parameters(LATITUDE)
    wave_params = calculate_wave_parameters(params, f)
    cr, ci = calculate_phase_speed(params, wave_params)
    
    # Calculate growth rate and e-folding time
    k = wave_params['k']
    nu = k * ci
    e_folding_time_days = np.log(2) / (nu * 3600 * 24)
    
    print(f'cr={cr:.4f} m/s  ||  ci={ci:.6f} m/s  ||  T={e_folding_time_days:.2f} days')
    
    # Get amplitudes and phases
    psi_amp, theta = calculate_streamfunction_amplitude(x3, params, wave_params, cr, ci)
    w_amp, sigma = calculate_vertical_velocity_amplitude(x3, params, wave_params, cr, ci, f)
    
    # Pre-calculate common terms
    l = wave_params['l']
    phase_x = k * (x1 - cr * t) + theta
    exp_growth = np.exp(nu * t)
    
    cos_x = np.cos(phase_x)
    sin_x = np.sin(phase_x)
    cos_y = np.cos(l * x2)
    sin_y = np.sin(l * x2)
    
    # Streamfunction
    psi = psi_amp * cos_x * exp_growth * sin_y
    
    # Geostrophic velocities
    ug = -exp_growth * l * psi_amp * cos_y * cos_x
    vg = -exp_growth * k * psi_amp * sin_y * sin_x
    
    # Ageostrophic velocities
    ua_term1 = l * cos_y * (-x2 * beta * cos_x + exp_growth * k**2 * psi_amp * sin_y)
    ua_term2 = k * sin_y * (cr * k * cos_x - nu * sin_x)
    ua = -(1 / f) * exp_growth * psi_amp * (ua_term1 + ua_term2)
    
    # Meridional ageostrophic velocity (expanded form for clarity)
    va_term1 = -2 * k * x2 * beta * sin_y * sin_x
    va_term2 = 2 * l * cos_y * (-nu * cos_x + cr * k * sin_x)
    va_term3 = exp_growth * k * l**2 * psi_amp * np.sin(2 * phase_x)
    va = (1 / (2 * f)) * exp_growth * psi_amp * (va_term1 + va_term2 + va_term3)
    
    # Horizontal divergence
    k2_plus_l2 = k**2 + l**2
    div_term1 = k2_plus_l2 * nu * cos_x * sin_y
    div_term2 = k * (cr * k2_plus_l2 + beta) * sin_x
    divergence = (1 / f) * exp_growth * psi_amp * sin_y * (div_term1 + div_term2)
    
    # Vorticity
    vort_term1 = exp_growth * k**2 * l**2 * psi_amp * (np.cos(2 * l * x2) - np.cos(2 * phase_x))
    vort_term2 = cos_x * (-l * beta * cos_y + k2_plus_l2 * (-f + x2 * beta) * sin_y)
    vorticity = (1 / f) * exp_growth * psi_amp * (vort_term1 + vort_term2)
    
    # Vertical velocity
    w = w_amp * np.cos(k * (x1 - cr * t) - sigma) * exp_growth * sin_y
    
    return psi, ug, vg, ua, va, w, divergence, vorticity
