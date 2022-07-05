# Structural-stability
This program performs the structural stability of a structure.

The user can perform the analysis of a structure and inspect the results.

The outflow of the program is as follows:

1. **Main.m** is the main routine. It calls all the necessary sub-routines and plot the results;

2. **Initial_calculations.m** does the calculations that only have to be done once, before the analysis starts;

3. **Stiffness_matrix.m** calculates the elementary stiffness matrix. This matrix is essential in the Newton-Raphson iterative scheme in order to achieve the solution of the problem.

4. **Deformed_config.m** plots the deformed configuration. It also allows the plot of the heat map of stresses/strains.

## Examples of the output of the program

* Load-displacement graphs:

![This is an image](Hat_Load-displacement.png)

![This is an image](LippedChannel_Load-displacement.png)

* Deformed configurations:

![This is an image](Hat_elastic.png)

![This is an image](LippedChannel_elastic.png)

![This is an image](LippedChannel_plastic.png)

* Heat map of stresses/strains:

![This is an image](LippedChannel_Strains.png) 
