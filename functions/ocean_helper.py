import numpy as np
import gsw
import matplotlib.pyplot as plt

def create_empty_ts(T_extent, S_extent, ax, p_ref = 0):
    '''
    Adapted from unisacsi Ocean.py (https://github.com/UNISvalbard/unisacsi)

    Creates an empty TS-diagram to plot data into.
    Parameters
    ----------
    T_extent : (2,) array_like
        The minimum and maximum conservative temperature.
    S_extent : (2,) array_like
        The minimum and maximum absolute salinity.
    p_ref : int, optional
        Which reference pressure to use. The following options exist:\n
        0:    0 dbar\n
        1: 1000 dbar\n
        2: 2000 dbar\n
        3: 3000 dbar\n
        4: 4000 dbar\n
        The default is 0.
    Returns
    -------
    None.
    '''
    sigma_functions = [gsw.sigma0,gsw.sigma1,gsw.sigma2,gsw.sigma3,gsw.sigma4]
    T = np.linspace(T_extent[0], T_extent[1], 100)
    S = np.linspace(S_extent[0], S_extent[1], 100)

    T,S = np.meshgrid(T,S)

    SIGMA = sigma_functions[p_ref](S,T)

    cs = ax.contour(S,T,SIGMA,colors='k',linestyles='--')
    ax.clabel(cs,fmt = '%1.2f')

    ax.set_ylabel('Conservative Temperature (°C)')
    ax.set_xlabel('Absolute Salinity')#[g kg$^{-1}$]
    #ax.set_title('$\Theta$ - $S_A$ Diagram')
    if p_ref > 0:
        ax.set_title('Density: $\sigma_{'+str(p_ref)+'}$',loc='left',fontsize=10)

    return


def classify_water_mass(theta, S, return_numeric=False):
    """
    Classify water masses into AW (Atlantic Water), PW (Polar Water),
    MW (Meltwater), or Mixed based on simple T/S thresholds.

    Parameters
    ----------
    theta : array-like
        Potential temperature [°C]
    S : array-like
        Salinity [psu]

    Returns
    -------
    labels : ndarray of str
        'AW', 'PW', 'MW', or 'Mixed'
    """

    theta = np.asarray(theta)
    S = np.asarray(S)

    labels = np.full(theta.shape, 'Mixed', dtype=object)

    # --- thresholds: adjust to your Fram Strait T–S space ---
    # Meltwater: very fresh regardless of temp
    mw_mask = S < 33.5
    labels[mw_mask] = 'MW'

    # Atlantic Water: warm & salty
    aw_mask = (S >= 34.8) & (theta >= 2.0)
    labels[aw_mask] = 'AW'

    # Polar Water: cold & relatively fresh
    pw_mask = (S <= 34.5) & (S>=33.5) & (theta <= 1.0)
    labels[pw_mask] = 'PW'

    # Resolve overlaps: MW overrides everything if very fresh
    labels[mw_mask] = 'MW'

    if return_numeric:
        label_to_code = {'Mixed': 0, 'AW': 1, 'PW': 2, 'MW': 3}
        codes = np.vectorize(label_to_code.get)(labels)
        return codes

    return labels