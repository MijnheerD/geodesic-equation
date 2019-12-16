# geodesic-equation

 This is the group project of Emilie Grégoire and Mitja Desmet for the course of computational physics in first master of physics at the VUB.

The structure is as follows:
- _project.py_ contains all definitions necessary to numerically integrate the geodesic equation and write the results to a file.
- _animation.py_ is used to plot the results from the integration using VPython. Best is to open this file in VIdle, the IDE which is installed alongside VPython. The functions are written specifically for the metrics we analysed to more accurately represent the situation, but the one for _Minkowski_ should be good for every metric.
- _testing.py_ was used to actually do all the numerical analysis and testing.
- _einstein.py_ uses the EinsteinPy package do also calculate geodesics. We used this to compare our results to.

All _.txt_ files are results from the different computations. The ones denote with _solveEinstein_ come from the EinsteinPy package. The ones which start with _solveGE_ come from calculations done with the solveGE function (see _project.py_). Similarly, the ones with _RK4_ are the ones which come from RK4 (also see _project.py_ for implementation).
- If there appears _Minskowski_ in the name, the considered spacetime was Minkowski spactime. Otherwise it was Schwarzschild spacetime. For this last one we always considered a mass of 5.972e24 kg. 
- The files with _phi_ or _theta_ are results from calculations with as initial condition only had a velocity in the phi or theta direction respectively.
- _mass_ refers to the only calculation with a different mass: the mass considerd was 5.972e30 kg.
- Finally _randphi_ had as initial condition a velocity in both the r and phi direction.

To end, _interesting-links_ contains some links to resources we used for our project.
