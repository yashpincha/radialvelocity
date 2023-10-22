import exoplanet as xo

periods = [] # FIX ME! include period of binaries
period_errs = [] # FIX ME! errors in the period

# semi amplitude = vmax - vmin / 2
Ks = xo.estimate_semi_amplitude(periods, x, y, yerr, t0s=None )
print(Ks, "m/s")
