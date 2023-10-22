import pymc3 as pm
import pymc3_ext as pmx
import aesara_theano_fallback.tensor as tt
import numpy as np

with pm.Model() as model:
    # Gaussian priors based on transit data
    transit_timings = pm.Normal("transit_timings", mu=np.array(t0s), sd=np.array(t0_errs), shape=2)
    log_period = pm.Normal(
        "log_period",
        mu=np.log(periods),
        sd=np.array(period_errs) / np.array(periods),
        shape=2,
        testval=np.log(periods),
    )
    period = pm.Deterministic("period", tt.exp(log_period))

    # Wide log-normal prior for semi-amplitude
    log_semi_amplitude = pm.Normal(
        "log_semi_amplitude", mu=np.log(Ks), sd=2.0, shape=2, testval=np.log(Ks)
    )

    # Eccentricity & argument of periasteron
    eccentricity_coords = pmx.UnitDisk("eccentricity_coords", shape=(2, 2), testval=0.01 * np.ones((2, 2)))
    eccentricity = pm.Deterministic("eccentricity", tt.sum(eccentricity_coords ** 2, axis=0))
    argument_of_periasteron = pm.Deterministic("argument_of_periasteron", tt.arctan2(eccentricity_coords[1], eccentricity_coords[0]))
    xo.eccentricity.vaneylen19("eccentricity_prior", multi=True, shape=2, fixed=True, observed=eccentricity)

    # Jitter & a quadratic RV trend
    noise = pm.Normal("noise", mu=np.log(np.median(yerr)), sd=5.0)
    radial_velocity_trend = pm.Normal("radial_velocity_trend", mu=0, sd=10.0 ** -np.arange(3)[::-1], shape=3)

    # Define the orbit
    orbital_parameters = xo.orbits.KeplerianOrbit(period=period, t0=transit_timings, ecc=eccentricity, omega=argument_of_periasteron)

    # Function for computing the full RV model
    def compute_radial_velocity_model(t, name=""):
        # First the RVs induced by the planets
        radial_velocity = orbital_parameters.get_radial_velocity(t, K=tt.exp(log_semi_amplitude))
        pm.Deterministic("radial_velocity" + name, radial_velocity)

        # Define the background model 
        # A Vandermonde matrix is a matrix with the terms of a geometric progression in each row. 
        A = np.vander(t - x_reference, 3)
        background = pm.Deterministic("background" + name, tt.dot(A, radial_velocity_trend))

        # Sum over planets and add the background to get the full model
        return pm.Deterministic("rv_model" + name, tt.sum(radial_velocity, axis=-1) + background)

    # Define the RVs at the observed times
    observed_rv_model = compute_radial_velocity_model(x)

    # Also define the model on a fine grid as computed above (for plotting)
    predicted_rv_model = compute_radial_velocity_model(t, name="pred")

    # Add the observation model
    error = tt.sqrt(yerr ** 2 + tt.exp(2 * noise))
    pm.Normal("observed_rv", mu=observed_rv_model, sd=error, observed=y)
