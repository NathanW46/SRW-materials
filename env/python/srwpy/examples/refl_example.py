try: #OC15112022
    import sys
    sys.path.append('../')
    from srwlib import *
    from uti_refl import *
except:
    from srwpy.srwlib import *
    from srwpy.uti_refl import *

_n_ph_en=100
_n_ang=100
_n_comp=1
_ph_en_start=5000
_ph_en_fin=30000
_ph_en_scale_type="log"
_ang_start=0.0001
_ang_fin=0.015
_ang_scale_type="lin"

Pt_refl=get_refl_arr('Pt', _n_ang=_n_ang, _n_ph_en=_n_ph_en, _n_comp=_n_comp, _ph_en_start=_ph_en_start, _ph_en_fin=_ph_en_fin, _ph_en_scale_type=_ph_en_scale_type,
                  _ang_start=_ang_start, _ang_fin=_ang_fin, _ang_scale_type=_ang_scale_type)

plot_refl_curves('Au', ph_en_start=1000, polarization='s', ph_en_fin=10000, ang_start=0.0001, ang_fin=17.4, suppressed=True)
plt.show()

plot_refl3d(
refl=Pt_refl,
_ang_start=_ang_start, _ang_fin=_ang_fin, _n_ang=_n_ang,
_ph_en_start=_ph_en_start, _ph_en_fin=_ph_en_fin, _n_ph_en=_n_ph_en,
n_comp=_n_comp,
)
a = get_refl_arr('Pt', _n_ang=2, _n_ph_en=2, _n_comp=1, _ph_en_start=_ph_en_start, _ang_start=_ph_en_fin)
print(mirror_reflectivity('Pt', energy=_ph_en_start, theta=_ang_start, polarization='s'))
print(a)
# plt.show()

                            