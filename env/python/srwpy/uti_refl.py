"""Reflectivity Utilities Module


Modules:

    calc_relf

.. moduleauthor:: Nathan Whittington
"""
# updated 10062025

from __future__ import print_function  # Python 2.7 compatibility
try:
    from xraydb.utils import (PLANCK_HC)
    from xraydb import (add_material, get_material, xray_delta_beta, mirror_reflectivity)
    xraydb_found = True
except:
    print("Warning: 'xraydb' not found. Defaulting Reflectivity to 1.")
    xraydb_found = False


import numpy as np
import matplotlib.pyplot as plt
from array import array
from scipy import constants




def get_refl_arr(
    _mat=None,
    _n_ph_en=1,
    _n_ang=1,
    _n_comp=1,
    _ph_en_start=0,
    _ph_en_fin=0,
    _ph_en_scale_type="lin",
    _ang_start=0,
    _ang_fin=0,
    _ang_scale_type="lin",
    _dens=None,
    _thickness=None,
):
    #TODO add docstring
    '''
    returns a C-aligned flat array complex array vs photon energy vs grazing angle vs component (sigma, pi) for the specified material(s)
    Example for 2 Photon Energy, 2 Angle, and 2 components:
         [(E_0, A_0, S), (E_1, A_0, S), (E_0, A_1, S), (E_1, A_1, S), (E_0, A_0, P), (E_1, A_0, P), (E_0, A_1, P), (E_1, A_1, P)]
    '''

    # Perfect reflectivity if xraydb is not installed
    if not xraydb_found:
        return 1
    

    nTot = int(_n_ph_en * _n_ang * _n_comp * 2)

    # photon energy sampling
    if _ph_en_scale_type == "lin":
        ph_en = np.linspace(_ph_en_start, _ph_en_fin, _n_ph_en)
    elif _ph_en_scale_type == "log":
        ph_en = np.logspace(np.log10(_ph_en_start), np.log10(_ph_en_fin), _n_ph_en)
    else:
        raise Exception(
            "Invalid photon energy scale type. Use 'lin' or 'log'.")

    # grazing angle sampling
    if _ang_scale_type == "lin":
        ang = np.linspace(_ang_start, _ang_fin, _n_ang)
    elif _ang_scale_type == "log":
        ang = np.logspace(np.log10(_ang_start), np.log10(_ang_fin), _n_ang)
    else:
        raise Exception("Invalid angle scale type. Use 'lin' or 'log'.")

    if isinstance(_mat, (list, np.ndarray, array)):
        refl = calc_multilayer_refl(
                        _mat, ang, ph_en, densities=_dens, thicknesses=_thickness, n_comp=_n_comp)
    else:
        refl = calc_refl(
                        _mat, ang, ph_en, density=_dens, n_comp=_n_comp)
    if nTot != len(refl):
        raise Exception(f"reflectivity array length {len(refl)} does not match expected length {nTot}.")
    return refl



def calc_refl(formula, ang, energy, density=None, n_comp=1, polarization='s'):
    """mirror reflectivity for a thick, single-layer mirror.

    Args:
       formula (string):            material name or formula ('Si', 'Rh', 'silicon')
       ang (float or nd-array):     mirror angle in radians
       energy (float or nd-array):  X-ray energy in eV
       density (float or None):     material density in g/cm^3
       n_comp (1 or 2):             number of polarization components (1 = s, 2 = s and p)
       polarization ('s' or 'p'):   mirror orientation relative to X-ray polarization (only works for single value reflectivity)
    Returns:
       real and imaginary mirror reflectivity values
    """
    # Edited from xraydb
            
    if density is None:
        if get_material(formula) == None:
            raise Exception(f"{formula} not found.\n" \
            "Specify Density or use uti_refl.add_mat() to add a custom material")
        formula, density = get_material(formula)

    delta, beta, _ = xray_delta_beta(formula, density, energy)
    n = 1 - delta - 1j*beta

    if isinstance(energy, (np.ndarray, list, array)) and isinstance(ang, (np.ndarray, list, array)):
        _ang, _energy = np.meshgrid(ang, energy, indexing='ij')
        _energy = _energy.ravel()
        _ang = _ang.ravel()
        n = np.tile(n, len(ang))

        # kiz is k in air/vacuum,  with n = 1.
        # ktz is k in mirror material, with n < 1.
        qf  = (2*np.pi/PLANCK_HC) * _energy
        kiz = qf * np.sin(_ang)
        ktz = qf * np.sqrt(n**2 - np.cos(_ang)**2)

        r = (kiz - ktz)/(kiz + ktz)
        r_flat = r.view(np.float64).ravel()

        if n_comp == 2:
            # p polarization
            ktz_p = ktz / n
            r_p = (kiz - ktz_p)/(kiz + ktz_p)
            
            r_flat_p = r_p.view(np.float64).ravel()
            r_flat = np.concatenate((r_flat, r_flat_p))
        elif n_comp != 1:
            raise Exception(f"n_comp must be 1 or 2, not {n_comp}.")

        print(f"len(r) = {r.shape}, n_comp = {n_comp}, n_ang = {len(ang)}, n_en = {len(energy)}")
        return array('d', r_flat)

    else:
        # if single values

        # kiz is k in air/vacuum,  with n = 1.
        # ktz is k in mirror material, with n < 1.
        qf  = 2*np.pi * energy/PLANCK_HC
        kiz = qf * np.sin(ang)
        ktz = qf * np.sqrt(n**2 - np.cos(ang)**2)
        if polarization == 'p':
            ktz = ktz / n
        r = (kiz - ktz)/(kiz + ktz)

        if isinstance(r, float):
            print("single value reflection")
            return (r*r.conjugate()).real
        else: 
            raise Exception(f"ang and energy must be nd-arrays, not {type(ang)} and {type(energy)}")


def calc_multilayer_refl(mats, ang, ph_en, thicknesses, densities=None, n_comp=1):
    """
    Calculate the reflectivity for multilayer surfaces.
    Args:
        mats (list): List of material names for each layer (top to bottom).
        ang: Incident angle in radians.
        ph_en: Photon energy in eV.
        thicknesses (list): List of layer thicknesses in Angstroms
        densities (list, optional): List of layer densities in g/cm^3.
        n_comp (int, optional): Number of polarization components (1 or 2).
    Returns:
        c-aligned array: Reflectivity values for the multilayer surface.
    """
    if thicknesses is None:
        raise Exception(f'Please provide thicknesses in Angstroms of each layer (exluding substrate)')
    if len(mats) != len(thicknesses)+1:
        raise Exception(f'number of materials ({len(mats)}) should match number of thicknesses excluding substrate ({len(thicknesses)})')
    if densities is not None and len(mats) != len(densities):
        raise Exception(f"If not None, number of densities ({len(densities)}) should match number of materials({len(mats)})")
    # if isinstance(ph_en, (list, np.ndarray, array)) or isinstance(ang, (list, np.ndarray, array)):
    #     raise Exception("Error: multiple energies or angles not supported yet")

    if densities is None:
        densities = [None]*len(mats)   
    refl = []

    k0 = 2 * np.pi * ph_en / PLANCK_HC
    n_layers = len(mats) 
    thicknesses.append(None) # substrate Layer

    kz = [k0*np.sin(ang)]  # air/vacuum layer (n = 1)
    kx = k0*np.cos(ang)
    d = [10]


    for i in range(n_layers):
        if densities[i] is None:
            if get_material(mats[i]) == None:
                raise Exception(f"{mats[i]} not found.\n" \
                "Specify Density or use uti_refl.add_mat() to add a custom material")
            mats[i], densities[i] = get_material(mats[i])

        delta, beta, _ = xray_delta_beta(mats[i], densities[i], ph_en)
        print(delta, beta)
        n_i = 1 - delta + 1j*beta
        # kz_i = k0*np.sqrt(n_i**2 - np.cos(ang)**2) 
        # print( kz_i)
        kz_i = np.sqrt((n_i*k0)**2 - kx**2 - 2*delta*k0**2 + 2j*beta*k0**2)
        print( kz_i)
        kz.append(kz_i)
        d.append(thicknesses[i])



    print("Layers in Order:")
    for i in range(len(kz)):
        print(f'Layer {i}, kz: {kz[i]}, d: {d[i]}')
    print('\n\n')
    r = (kz[-2] - kz[-1])/(kz[-2] + kz[-1])
    print(f'refl of substrate = {(r*r.conjugate()).real}')
    # exclude substrate >----------vv
    for i in reversed(range(len(kz)-2)):
        print(f"Layer {i}: {kz[i]}, {d[i]}")
        fresnel_r = (kz[i] - kz[i+1])/(kz[i] + kz[i+1])
        print("r' =", fresnel_r)
        p = np.exp(2j*d[i+1]*kz[i+1])
        print('Phase shift =', p)
        r = (fresnel_r + r*(p**2))/(1 + fresnel_r*r*(p**2))
        print(f"refl between Layer {('air' if i == 0 else i)} and {'substrate' if (i+2) == len(kz) else (i+1)} is {(r*r.conjugate()).real}")


    # refl = np.empty(r.size*n_comp*2, dtype=np.float64)
    # if n_comp == 1:
    #     refl[0::2] = r.real
    #     refl[1::2] = r.imag

    # print(f'final refl of {mats} mirror = {refl}')
    # return array('d', refl)
    return (r*r.conjugate()).real


def add_mat(name: str, formula: str, density: float, categories: list[str] | None = None) -> None:
    """add a material to the users local material database

    Args:
        name (str): name of material
        formula (str): chemical formula
        density (float): density
        categories (list of strings or None): list of category names

    Returns:
        None

    Notes:
        the data will be saved to the file 'xraydb/materials.dat' in the
        configuration folder, and will be useful in subsequent sessions.


    Examples:
        >>> xraydb.add_material('boron carbide', 'B4C', 2.5, categories=['ceramic'])
    """

    if density == None:
        raise Exception("Error: density of {formula} not found")
    else:
        add_material(name, formula, density, categories)
    return
    

def plot_refl_curves(mat, **kwargs):
    '''
    - mat: material (provide density if the material is not defined)
    - kwargs: 
        - density (g/cm^3)
        - ang_start = 0.1 mrad
        - ang_fin = 15 mrad
        - ph_en_start = 1000 eV
        - ph_en_fin = 35000 eV
        - polarization = 's'
    '''

    try:
        from matplotlib.widgets import Slider
    except:
        print("Error: matplotlib.widgets.Slider not found")
        return

    #TODO check for material
    density = kwargs.get('density', None)
    ang_start = kwargs.get('ang_start', 0.1)
    ang_fin = kwargs.get('ang_fin', 15)
    ph_en_start = kwargs.get('ph_en_start', 1000)
    ph_en_fin = kwargs.get('ph_en_fin', 35000)
    polarization = kwargs.get('polarization', 's')
    suppressed = kwargs.get('suppressed', False)

    # initial params
    init_ang = 2
    init_ph_en = 10000

    ph_en_range = np.linspace(ph_en_start, ph_en_fin, 200)
    ang_range = np.linspace(ang_start, ang_fin, 50)

    fig = plt.figure(figsize=(10, 4.75))

    # Left plot and slider
    ax1 = fig.add_axes([0.07, 0.15, 0.375, 0.7])
    ax_ang_slider = fig.add_axes([0.475, 0.2, 0.02, 0.5])

    # Right plot and slider
    ax2 = fig.add_axes([0.565, 0.15, 0.375, 0.7])
    ax_ph_en_slider = fig.add_axes([0.515, 0.2, 0.02, 0.5])

    # initial plot
    if polarization == 's':
        ph_en_line, = ax1.plot(ph_en_range, mirror_reflectivity(mat, init_ang/1000, ph_en_range, density))
        ang_line, = ax2.plot(ang_range, mirror_reflectivity(mat, ang_range/1000, init_ph_en, density))
    elif polarization == 'p':
        ph_en_line, = ax1.plot(ph_en_range, mirror_reflectivity(mat, init_ang/1000, ph_en_range, density, polarization='p'))
        ang_line, = ax2.plot(ang_range, mirror_reflectivity(mat, ang_range/1000, init_ph_en, density, polarization='p'))
    else:
        raise Exception(f"Invalid polarization '{polarization}'. Use 's' or 'p'.")

    fig.suptitle(f'{mat} Reflectivity Curves: {polarization} polarization', fontsize=18)

    #slider text
    fig.text(0.47, 0.5, 'Grazing Angle (mrad)', va='center', ha='right', fontsize=10, rotation=90)
    fig.text(0.555, 0.5, 'Photon Energy (eV)', va='center', ha='right', fontsize=10, rotation=90)

    ax1.set_title('Reflectivity vs. Photon Energy', fontsize=15)
    ax1.set_xlabel('Photon Energy (eV)')
    ax1.set_ylabel('Reflectivity')

    ax2.set_title('Reflectivity vs. Grazing Angle', fontsize=15)
    ax2.set_xlabel('Grazing Angle (mrad)')
    ax2.set_ylabel('Reflectivity')
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()

    if not suppressed:
        ax1.set_ylim(0,1.05)
        ax2.set_ylim(0,1.05)

    # sliders
    ang_slider = Slider(ax_ang_slider, '', 0.5, 20, valinit=init_ang, orientation='vertical', valstep=0.5, )
    ph_en_slider = Slider(ax_ph_en_slider, '', 500, 50000, valinit=init_ph_en, orientation='vertical', valstep=50)

    def update(ang, ph_en):
        ph_en_line.set_ydata(mirror_reflectivity(mat, ang/1000, ph_en_range, density, polarization=polarization))
        ang_line.set_ydata(mirror_reflectivity(mat, ang_range/1000, ph_en, density, polarization=polarization))
        for ax in (ax1, ax2):
            ax.relim()
            ax.autoscale_view('y')
        fig.canvas.draw_idle()

    def on_timer(event):
        ang = ang_slider.val
        ph_en = ph_en_slider.val
        update(ang, ph_en)

    # Start timer
    timer = fig.canvas.new_timer(interval=250)
    timer.add_callback(on_timer, None)
    timer.start()


def plot_refl2d(
    refl: np.ndarray,
    _ang_start: float, _ang_fin: float, _n_ang: int,
    _ph_en_start: float, _ph_en_fin: float, _n_ph_en: int,
    n_comp: int = 1,
    constant: float = 10000,
    material: str = ""
):
    """
    Plot 2D reflectivity curves

    Parameters
    ----------
    refl : np.ndarray
        Complex C-aligned array of reflectivity values (photon energy vs grazing angle vs polarization)
    _ang_start : float
        Initial grazing angle in rad
    _ang_fin : float
        Final grazing angle in rad (exlusive)
    _n_ang : int
        Number of grazing angles
   _ph_en_start : float
        Initial photon energy in eV
    _ph_en_fin : float
        Final photon energy in eV (exlusive)
    _n_ph_en : int
        Number of photon energies
    n_comp : int, optional
        Number of polarization components. Defaults to 1.
    p_index : int, optional
        Index of the polarization component to plot (0 = s, 1 = p). Defaults to 0
    constant : float
        Constant value to plot against (either energy (eV) or angle (rad))
    """
    raise Exception("Not implemented yet")
    nTot = int(_n_ph_en * _n_ang * n_comp * 2)
    y = 0
    plt.figure(figsize=(8, 6))

    if constant > np.pi/2:
        x = np.linspace(_ang_start, _ang_fin, _n_ang)
        eStep = (_ph_en_fin - _ph_en_start)/(_n_ph_en) if _n_ph_en > 1 else 1
        e_i = round((constant - _ph_en_start)/eStep + 0.00001)
        id_start = 0 # (((e_i) * n_comp + p_index) * 2)
        id_fin =  _n_ang * 2 #(((e_i + _n_ph_en * x[-1]) * n_comp + p_index) * 2)
        print(id_start, id_fin)
        y = [((refl[i] + 1j*refl[i+1])*(refl[i] - 1j*refl[i+1])).real for i in range(id_start, id_fin, n_comp*2)]

        title = f"{material + ' ' if material != '' else ''}Reflectivity vs. Grazing Angle at {constant:.2f} eV"
        plt.xlabel("Ang (mrad)")

        plt.plot(x*1000, y)

    else:
        raise Exception("Not implemented")
        energy = np.linspace(_ph_en_start, _ph_en_fin, _n_ph_en)

    plt.ylabel("Reflectivity")
    plt.title(title)
    plt.tight_layout()

    

    



def plot_refl3d(
    refl: np.ndarray,
    _ang_start: float, _ang_fin: float, _n_ang: int,
    _ph_en_start: float, _ph_en_fin: float, _n_ph_en: int,
    n_comp: int = 1,
    p_index: int = 0,
    title: str = "Reflectivity"
):
    """
    Plot 3D surface of the real part of reflectivity for a given polarization component.

    Parameters:
    - refl: 1D float array Complex C-Style array vs energy vs angle vs polarization (s/p)
    - _ang_start, _ang_fin, _n_ang: angle range and resolution
    - _ph_en_start, _ph_en_fin, _n_ph_en: energy range and resolution
    - n_comp: number of polarization components (default 1 = s)
    - p_index: polarization index to plot (0 = s, 1 = p)
    - title: plot title
    """
    if not isinstance(refl, (np.ndarray, list, array)):
        raise Exception("Reflectivity maps require C-aligned array of reflectivity data")


    ang = np.linspace(_ang_start, _ang_fin, _n_ang)
    energy = np.linspace(_ph_en_start, _ph_en_fin, _n_ph_en)

    ang *= 1000

    T, E = np.meshgrid(ang, energy, indexing='ij')
    Z = np.zeros((_n_ang, _n_ph_en))

    for j in range(_n_ang):   
        for i in range(_n_ph_en):
            idx = (((i + _n_ph_en * j) * n_comp + p_index) * 2)
            Z[j, i] = refl[idx]

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(E, T, Z, cmap='viridis', edgecolor='none')

    ax.set_xlabel("Photon Energy (eV)")
    ax.set_ylabel("Angle (mrad)")
    ax.set_zlabel("Reflectivity")
    ax.set_title(title)

    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label="Reflectivity")
    plt.tight_layout()


# def get_indices_from_arr(
#     index,
#     _n_ph_en,
#     _n_ang,
#     _n_comp,
#     _ph_en_start,
#     _ph_en_fin,
#     _ang_start,
#     _ang_fin,
# ):
#     ang = np.linspace(_ang_start, _ang_fin, _n_ang)
#     energy = np.linspace(_ph_en_start, _ph_en_fin, _n_ph_en)
#     nTot = int(_n_ph_en * _n_ang * _n_comp * 2)

#     e_i = int(index//(nTot/_n_ph_en))
#     a_i = int(index%(nTot/_n_ang))//(2*_n_comp)
#     if _n_comp == 1:
#         p_i = 0
#     elif _n_comp == 2:
#         p_i = index%4 // 2
#     c_i = index%2

#     eStep = (_ph_en_fin - _ph_en_start)/(_n_ph_en - 1)
#     angStep = (_ang_fin - _ang_start)/(_n_ang - 1)
#     print(eStep, angStep)
#     print(energy[1] - energy[0], ang[1] - ang[0])

#     print(f'For index {index}: Energy[{e_i}] = {energy[e_i]} eV, Angle[{a_i}] = {ang[a_i]} mrad, Polarization = {p_i}, Component = {c_i}', )