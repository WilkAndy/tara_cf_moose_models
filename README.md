# Geomechanical response due to nonisothermal fluid injection into a reservoir

The repository contains MOOSE models to compare with and extend the LaForce et al. semi-analytic solutions of nonisothermal fluid injection into a reservoir.

The MOOSE models and results are described in:

- Green, C.P. and Wilkins, A. and Ennis-King, J. and LaForce, T., "Geomechanical response due to nonisothermal fluid injection into a reservoir", Advances in Water Resources (2021), https://doi.org/10.1016/j.advwatres.2021.103942

The LaForce et al. model and solutions are defined in:

- LaForce, T. and Ennis-King, J. and Paterson, L., "Semi-analytical solutions for nonisothermal fluid injection including heat loss from the reservoir: Part 1. Saturation and temperature", Advances in Water Resources 73 (2014) 227--234, https://www.doi.org/10.1016/j.advwatres.2014.08.008
- LaForce, T. and Mijic, A. and Ennis-King, J. and Paterson, L., "Semi-analytical solutions for nonisothermal fluid injection including heat loss from the reservoir: Part 2.  Pressure and stress", Advances in Water Resources 73 (2014) 242--253, https://www.doi.org/10.1016/j.advwatres.2014.08.009

The MOOSE models contained herein run with the following version:

```
Framework Information:
MOOSE Version:           git commit e900438bb7 on 2021-05-10
LibMesh Version:         27141d18f3137f77e33cdb3d565fd38ebfbfc46f
PETSc Version:           3.14.6
Current Time:            Mon May 17 08:24:23 2021
Executable Timestamp:    Mon May 17 08:05:13 2021
```

Download, install and usage instructions for MOOSE may be found at https://mooseframework.inl.gov/ .   The "benchmark" MOOSE model that is guaranteed to run with the latest version of MOOSE may be found at https://mooseframework.inl.gov/modules/porous_flow/thm_example.html (this model does not contain any of the extensions described in the Green et al. paper).

Finally, note two things:

- You will probably receive slightly different results from your own install of MOOSE compared with those in this repository.  This is probably due to you having a slightly different PETSc configuration than was used here.  The differences won't result in any noticable changes in the plots.
- Before running the python scripts to generate plots based on your results, you will probably have to modify those python scripts to ensure the correct CSV files are used, and the correct columns in those CSV files are plotted.



