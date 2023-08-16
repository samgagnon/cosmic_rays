import h5py
import numpy as np
from astropy import constants as const2
from astropy import units as u

pc = const2.pc.cgs.value
kpc = const2.kpc.cgs.value
Msol = const2.M_sun.cgs.value
Myr = ((1.*u.Myr).to(u.s).value)

def read_units(fname):
    print("reading file: ", fname)
    f = h5py.File(fname, "r")
    # get momentum boundaries from file
    try:
        UnitMass = f[u'Parameters'].attrs.get("UnitMass_in_g")
        UnitLen  = f[u'Parameters'].attrs.get("UnitLength_in_cm")
        UnitVelo = f[u'Parameters'].attrs.get("UnitVelocity_in_cm_per_s")
        BoxSize  = f[u'Parameters'].attrs.get("BoxSize")
    except:
        print("Unit parameters not found!!!")
    f.close()
    return UnitMass, UnitLen, UnitVelo, BoxSize

def compute_arepo_units(UnitMass, UnitLength, UnitVelocity, print_units=False):
    UnitTime = UnitLength/UnitVelocity
    d = { "UnitLength"   : UnitLength,
          "UnitMass"     : UnitMass,
          "UnitVelocity" : UnitVelocity,
          "UnitTime"     : UnitLength/UnitVelocity,
          "UnitDensity"  : UnitMass/(UnitLength**3),
          "UnitColdens"  : UnitMass/(UnitLength**2),
          "UnitEnergy"   : UnitMass*UnitVelocity*UnitVelocity,
          "UnitEspecific": UnitVelocity*UnitVelocity,
          "UnitEdensity" : UnitMass / UnitLength / UnitTime**2,
          "UnitPressure" : UnitMass / UnitLength / UnitTime**2,
          "UnitBfield"   : np.sqrt(UnitMass) / np.sqrt(UnitLength) / UnitTime,
          "UnitGrav"     : UnitLength**3 / (UnitMass * UnitTime**2),
          "UnitPotential": UnitLength**2 / (UnitTime**2),
          # non cgs
          "UnitColdens_Msol_pc2" : UnitMass/(UnitLength**2) / Msol * pc*pc

    }
    if print_units:
        print()
        print("Arepo units")
        # find longest key name
        maxlen = len(max(d, key=len))
        for key in d:
            print(key.ljust(maxlen), ":", d[key])
        print()
    return d


def read_momentum_bins(fname):
    print("reading file: ", fname)
    f = h5py.File(fname, "r")
    Nspec = f[u'PartType0/CRspecEnergy'].shape[1]
    # get momentum boundaries from file
    try:
        pmin = f[u'Parameters'].attrs.get("CRspec_pmin")
        p0   = f[u'Parameters'].attrs.get("CRspec_p0")
        p1   = f[u'Parameters'].attrs.get("CRspec_p1")
        pmax = f[u'Parameters'].attrs.get("CRspec_pmax")

        pf = np.zeros(Nspec+1)
        pf[0]  = pmin
        pf[-1] = pmax
        pf[1:-1] = np.logspace(np.log10(p0), np.log10(p1), Nspec-1, endpoint=True)

        pi   = np.sqrt(pf[:-1]*pf[1:])
    except:
        print("momentum parameters not found!!!")
    return pf, pi