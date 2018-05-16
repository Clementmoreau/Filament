# Filament

This set of files is provided as supplementary material to the article "The asymptotic coarse-graining formulation of slender-rods, bio-filaments and flagella", to be published in Journal of the Royal Society Interface, authored by C. Moreau, L. Giraldi, and H. Gadêlha.
arxiv link : https://arxiv.org/abs/1712.02697

These files can be run on MATLAB in order to simulate the behaviour of an elastic filament in a fluid at low Reynolds number.

How to use it :

MAIN file : choose 
- final time and time step
- "sperm number" (the lower, the more rigid the filament is)
- number of segments (the simulation is more accurate with more segments but don't use more than approx. 80 - a good choice is between 20 and 40)
- an initial shape for the filament
- relaxation, oscillation or magnetic case

The code solves the dynamics and returns the time (time) and state (traj) vectors.
A graphic visualisation starts at the end. 

If any questions contact me at clement dot moreau at inria dot fr

