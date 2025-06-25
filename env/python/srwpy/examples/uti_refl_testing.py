try: #OC15112022
    import sys
    sys.path.append('../')
    from srwlib import *
    from uti_plot import *
    from uti_refl import *
except:
    from srwpy.srwlib import *
    from srwpy.uti_plot import *
    from srwpy.uti_refl import *
import os
import time
import copy
import numpy as np
import xraydb

ph_en = 10000
ang = 0.005
refl=calc_multilayer_refl(['Au', 'Si'], ph_en=ph_en, ang=ang, thicknesses=[10])


print(f"Single mirror refl of Au at {ph_en/1000} keV and {ang*1000} is {mirror_reflectivity('Au', ang, ph_en)}")
print(f"Single mirror refl of Si at {ph_en/1000} keV and {ang*1000} is {mirror_reflectivity('Si', ang, ph_en)}")


# print('Testing uti_materials module')

# def mirror_reflectivity(formula, theta, energy, density=None,
#                         roughness=0.0, polarization='s'):
#     """mirror reflectivity for a thick, singl-layer mirror.

#     Args:
#        formula (string):           material name or formula ('Si', 'Rh', 'silicon')
#        theta (float or nd-array):  mirror angle in radians
#        energy (float or nd-array): X-ray energy in eV
#        density (float or None):    material density in g/cm^3
#        roughness (float):          mirror roughness in Angstroms
#        polarization ('s' or 'p'):  mirror orientation relative to X-ray polarization

#     Returns:
#        mirror reflectivity values

#     Notes:
#        1. only one of theta or energy can be an nd-array
#        2. density can be `None` for known materials
#        3. polarization of 's' puts the X-ray polarization along the mirror
#           surface, 'p' puts it normal to the mirror surface. For
#           horizontally polarized X-ray beams from storage rings, 's' will
#           usually mean 'vertically deflecting' and 'p' will usually mean
#           'horizontally deflecting'.
#     """
#     from xraydb.materials import get_material
#     if density is None:
#         formula, density = get_material(formula)

#     delta, beta, _ = xray_delta_beta(formula, density, energy)
#     n = 1 - delta - 1j*beta

#     # kiz is k in air/vacuum,  with n = 1.
#     # ktz is k in mirror material, with n < 1.
#     qf  = 2*np.pi * energy/PLANCK_HC
#     kiz = qf * np.sin(theta)
#     ktz = qf * np.sqrt(n**2 - np.cos(theta)**2)

#     # polarization correction will be tiny for small angles
#     if polarization == 'p':
#         ktz = ktz / n

#     r_amp = (kiz - ktz)/(kiz + ktz)
#     if roughness > 1.e-12:
#         r_amp = r_amp * np.exp(-2*(roughness**2*kiz*ktz))
#     return (r_amp*r_amp.conjugate()).real


# #******************* Reflectivity Parameters

# _n_ph_en=100
# _n_ang=100
# _n_comp=2
# _ph_en_start=1000
# _ph_en_fin=10000
# _ph_en_scale_type="lin"
# _ang_start=0.01
# _ang_fin=0.1
# _ang_scale_type="lin"

# nTot = int(_n_ph_en * _n_ang * _n_comp * 2)
# print(f"nTot = {nTot} = {_n_ph_en} * {_n_ang} * {_n_comp} * 2")


# eStep = (_ph_en_fin - _ph_en_start)/(_n_ph_en-1)
# angStep = (_ang_fin - _ang_start)/(_n_ang-1)

# polar = {0:'s', 1:'p'}

# test_ang = _ang_start + angStep*99
# test_en= _ph_en_start+ eStep*99
# test_polar = 0

# ie = round((test_en - _ph_en_start)/eStep + 0.00001)
# iang = round((test_ang - _ang_start)/angStep + 0.00001)
# ipol = test_polar
# i = (2*(_n_comp*(ie*_n_ang + iang)+ipol))
# print(f"ie={ie}, iang={iang}, ipol={ipol}, i={i}")
# formula, density = xraydb.get_material('B4C')
# delta, beta, _ = xray_delta_beta(formula, density, energy=test_en)
# n = 1 - delta - 1j*beta

# qf  = (2*np.pi/PLANCK_HC) * test_en
# kiz = qf * np.sin(test_ang)
# ktz = qf * np.sqrt(n**2 - np.cos(test_ang)**2)

# r_amp = ((kiz - ktz)/(kiz + ktz))
# r_amp = np.abs(r_amp)**2


# refl=get_refl_arr(formula, _n_ang=_n_ang, _n_ph_en=_n_ph_en, _n_comp=_n_comp, _ph_en_start=_ph_en_start, _ph_en_fin=_ph_en_fin, _ph_en_scale_type=_ph_en_scale_type,
#                   _ang_start=_ang_start, _ang_fin=_ang_fin, _ang_scale_type=_ang_scale_type, _dens=density)

# print(f"refl from array refl[{i}] =", np.abs(refl[i] + 1j*refl[i+1])**2)
# recalc_refl = mirror_reflectivity(formula, energy=test_en, theta=test_ang, polarization=polar[test_polar], density=density)

# print(f"recalculated refl for {formula} at the same index is", recalc_refl)

# print("checking every 100th index")
# mismatch = False
# i = 0
# # for i in range(0, 405):
# for i in range(0, len(refl)-1, 100):
#     e_i = int(i//(nTot/_n_ph_en))
#     a_i = int(i%(nTot/_n_ang))//(2*_n_comp)
#     if _n_comp == 1:
#         p_i = 0
#     elif _n_comp == 2:
#         p_i = i%4 // 2
#     c_i = i%2
#     if c_i == 1:
#         #if imaginary
#         continue
#     recalc_refl_check = mirror_reflectivity(formula, energy=_ph_en_start+eStep*e_i, theta=_ang_start+angStep*a_i, polarization=polar[p_i])
#     if not np.isclose(np.abs(refl[i] + 1j*refl[i+1])**2, recalc_refl_check, rtol=0.001):
#         print(f"i={i}, e_i={e_i}, a_i={a_i}, p_i={p_i}, c_i={c_i}", end=' ')
#         print(np.abs(refl[i] + 1j*refl[i+1])**2, recalc_refl_check)
#         mismatch = True
# if not mismatch:
#     print('All tested values match within 0.1% tolerance')

# i = 0
# print("checking every 48th index using for loops")
# for i in range(0, _n_ph_en, 4):
#     for j in range(0, _n_ang, 3):
#         for k in range(0, _n_comp):
#             recalc_refl_check = mirror_reflectivity(formula, energy=_ph_en_start+eStep*i, theta=_ang_start+angStep*j, polarization=polar[k])
#             index = ((i*_n_comp*_n_ang) + (j*_n_comp) + k)*2
#             if not np.isclose(np.abs(refl[index] + 1j*refl[index+1])**2, recalc_refl_check, rtol=0.001):
#                 print(f"i={index}, e_i={i}, a_i={j}, p_i={k}, c_i={1}", end=' ')
#                 print(np.abs(refl[index] + 1j*refl[index+1])**2, recalc_refl_check)
#                 mismatch = True
# if not mismatch:
#     print('All tested values match within 0.1% tolerance')


# print('All tested values match within 0.1% tolerance using for loops')
