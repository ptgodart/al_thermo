Aluminum-Water Reaction Thermodynamics
======================================

Code for computing Gibbs free energy for different Al reaction product species.
The raw enthalpy values for the different species at different temperatures are
tabulated in the data folder. One file for each species. These values were taken
from experimental results published in literature.

How to run:
-----------

  1. Run `delta_G_using_numerical_Integration_of_H.m` to generate `delta_G`
     values for the different species at different temperatures.
  2. Run `delta_g_computation.m`, which takes `delta_G` values from previous
     step and applies stoichiometric ratios and the effects of pressure to
     compute final `delta_G` values for each reaction. This script also plots
     the values against temperature to illustrate where favorability shifts.
