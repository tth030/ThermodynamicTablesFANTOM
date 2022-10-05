
import numpy as np
from ..tools import *
from scipy.interpolate import RegularGridInterpolator

g = 9.80665 # (m/s2)

class strength_profile(object):
    
    def __init__(self,P=[],Z=[],stress_failure=[],stress_failure_weak=[],strain_rate=1e-16):
        self.stress_failure       = stress_failure     # pressure (Pa)
        self.stress_failure_weak  = stress_failure_weak     # pressure (Pa)
        self.Z  = Z     # depth (m)
        self.strain_rate = strain_rate # percentage
        self.dz = 200    # output sampling (m)

def compute_strength_profile(layers,geotherm,P_rho_DepthProfile,strain_rate=1e-16,printOut=True):
    """ Compute the strength profiles based on properties
    given by the list of layers for given temperature and pressure profile.
    Each "layer" is an object defined in core.py of geodyn1d """

    StrengthProfile = strength_profile()
    nlayers         = len(layers)
    StrengthProfile.strain_rate = strain_rate

    # 1) Browsing P,T values to compute the strength profile
    stress_failure      = np.zeros_like(geotherm.T).tolist()
    stress_failure_weak = np.zeros_like(geotherm.T).tolist()
    P                   = P_rho_DepthProfile.P.tolist()
    T                   = geotherm.T.tolist()
    Z                   = geotherm.Z.tolist()

    nl  = len(Z)
    #print('min {} max {} len {}'.format(min(Z),max(Z),len(Z)))
    #print('min {} max {} len {}'.format(min(T),max(T),len(T)))

    for ni in range(0,nl):
        z = Z[nl-ni-1]
        t = T[nl-ni-1]
        cumdepth = 0.0e0
        for ilayer in range(0,nlayers):
            #cumdepth = cumdepth + layers[ilayer].thickness
            if (z<cumdepth+ layers[ilayer].thickness):
                break
            cumdepth = cumdepth + layers[ilayer].thickness

        #print('{} {} {} {}'.format(z,cumdepth,layers[ilayer].thickness,ilayer))
        #layers[ilayer].display()

        cohesion      = layers[ilayer].material.mechanicalProperties.plasticity.cohesion1
        cohesion_weak = layers[ilayer].material.mechanicalProperties.plasticity.cohesion2
        phi           = layers[ilayer].material.mechanicalProperties.plasticity.internalAngleOfFriction1
        phi_weak      = layers[ilayer].material.mechanicalProperties.plasticity.internalAngleOfFriction2
        PLB           = layers[ilayer].material.mechanicalProperties.plasticity.PowerLawBreakdown
        f_power_law   = layers[ilayer].material.mechanicalProperties.flowlaw.Fc
        A_power_law   = layers[ilayer].material.mechanicalProperties.flowlaw.A
        n_power_law   = layers[ilayer].material.mechanicalProperties.flowlaw.n
        Ea_power_law  = layers[ilayer].material.mechanicalProperties.flowlaw.Ea
        Va_power_law  = layers[ilayer].material.mechanicalProperties.flowlaw.Va
        
        pressure    = P[nl-ni-1]
        temperature = T[nl-ni-1]

        if (ni==0):
            dz       = 0.0e0
            stress_failure[nl-ni-1]      = cohesion
            stress_failure_weak[nl-ni-1] = cohesion_weak
        else:
            dz= (Z[nl-ni-1] - Z[nl-ni])
            sigmaP      = pressure * np.sin(np.deg2rad(phi)) + cohesion * np.cos(np.deg2rad(phi))
            sigmaP_weak = pressure * np.sin(np.deg2rad(phi_weak)) + cohesion_weak * np.cos(np.deg2rad(phi_weak))

            if PLB :
                if sigmaP > PLB:
                    sigmaP = PLB

            sigmaV = f_power_law * A_power_law**(-1/n_power_law) * strain_rate**(1/n_power_law) * \
                     np.exp((Ea_power_law + pressure*Va_power_law)/(n_power_law*8.3145*temperature))

            stress_failure[nl-ni-1]      = min(sigmaP,sigmaV)
            stress_failure_weak[nl-ni-1] = min(sigmaP_weak,sigmaV)
        #print('{} {}'.format(Z[nl-ni-1],rho[nl-ni-1]))
       
    StrengthProfile.stress_failure       = np.asarray(stress_failure)
    StrengthProfile.stress_failure_weak  = np.asarray(stress_failure_weak)
    StrengthProfile.Z                    = np.asarray(Z)
    StrengthProfile.dz                   = dz

    return StrengthProfile

