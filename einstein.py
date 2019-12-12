import numpy as np
from astropy import units as u
from einsteinpy.coordinates import SphericalDifferential, CartesianDifferential
from einsteinpy.metric import Schwarzschild

M = 5.972e24 * u.kg
sph_coord = SphericalDifferential(306.0 * u.m, np.pi/2 * u.rad, -np.pi/6*u.rad,
                          0*u.m/u.s, 0*u.rad/u.s, 1900*u.rad/u.s)
obj = Schwarzschild.from_coords(sph_coord, M , 0* u.s)

end_tau = 0.01 # approximately equal to coordinate time
stepsize = 0.3e-6
ans = obj.calculate_trajectory(end_lambda=end_tau, OdeMethodKwargs={"stepsize":stepsize})
print(ans)
