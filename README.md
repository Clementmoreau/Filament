# Filament

These files can be run on MATLAB in order to simulate the behaviour of an elastic filament in a fluid at low Reynolds number.

How to use it :

MAIN file : choose
- final time
- "sperm number" (the lower, the more rigid the filament is)
- number of segments (the simulation is more accurate with more segments but don't use more than approx. 80 - a good choice is between 20 and 40)
- an initial shape for the filament
- if using second_member_mag, choose a time-varying magnetic field

The code solves the dynamics and returns a few plots to visualise it.

You can change the name of the solving function (line 54) to study different phenomenons:
- second_member_Nparam : standard free ends evolution
- second_member_Nparam_counterbend : adding a torque to model the counterbend phenomenon (due to a sliding displacement between two filaments at the base)
- second_member_Nparam_oscillation : forced angular actuation at the proximal end (the angle of actuation can be freely chosen)
- second_member_mag : adding an external magnetic field

Other MAIN files :

MAIN_buckling : buckling phenomenon with different initial conditions
MAIN_bundle : simulation for a double filament with springs linking them

If any questions contact me at clement dot moreau at inria dot fr

