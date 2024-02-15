from obspy.core import Stream
from obspy.signal.rotate import rotate_ne_rt
from telewavesim import utils as ut
from telewavesim import wiggle as wg
import numpy as np 

# LAYmodel,ID,ph,samprate,inc,synthperiod,nsamps,cutf,sourc
def run_telewavesim(modfile, wvtype, npts, dt, dp, use_obs, c, rhof, slow, baz): 
    # Back azimuth is in degrees: used in rotate_ne_rt. 
    model = ut.read_model(modfile)
    trxyz = ut.run_plane(model, slow, npts, dt, baz=baz, wvtype=wvtype,
                        obs=use_obs, dp=dp, c=c, rhof=rhof)

    ntr = trxyz[0] # North component
    etr = trxyz[1] # East component
    ztr = trxyz[2] # Vertical component

    # Copy to radial and transverse. Temporary: Will do rotation. 
    rtr = ntr.copy() # Radial component
    ttr = etr.copy() # Transverse component

    # Rotate to radial and transverse
    if not ((baz % 360) == 0):
        rtr.data, ttr.data = rotate_ne_rt(ntr.data, etr.data, baz)

    traces = np.array([rtr, ttr, ztr]) # TODO double check. do we want radial transverse? 
    tt = np.arange(0, npts)*dt # TODO double check. 
    # status = ''
    # cmdout = '' 
    return traces, tt # , status, cmdout


# Test the function
if __name__ == '__main__': 
    import matplotlib.pyplot as plt 
    modfile = './demo.txt'
    wvtype = 'P'

    npts = 3000 # Number of samples
    dt = 0.01   # Sample distance in seconds

    dp = 0 # 2000. # Deployment depth below sea level in meters
    # use_obs = dp > 0 # Only use OBS code if deployment depth is beneath sea level. 
    use_obs = False

    c = 1.5      # P-wave velocity in salt water (km/s)
    rhof = 1027. # Density of salt water (kg/m^3)

    slow = 0.06 # Horizontal slowness (or ray parameter) in s/km 
    baz = 0.    # Back-azimuth direction in degrees (has no influence if model is isotropic)

    traces, tt = run_telewavesim(modfile, wvtype, npts, dt, dp, use_obs, c, rhof, slow, baz)

    fig = plt.figure() 
    plt.ion()
    plt.plot(tt, traces.T) 
    plt.savefig('./figs_test/test1.jpeg') 