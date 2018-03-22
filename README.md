# Single-cell-imaging; software associated with paper "Induction of CD4 T cell memory by local cellular collectivity", Michal Polonsky, Jacob Rimer, Amos Kern-Perets, Irina Zaretsky, Stav Miller, Chamutal Bornstein, Eyal David, Naama Kopelman, Ziv Porat, Benny Chain, Nir Friedman 

The repository contains a zip file SingleCell_Sim.zip which contains a Matlab code file (.m) and an associated .mat file  Both files should be placed in the same folder. 

The code is a Matlab function that reads the attached .m file which contains the distributions that the function uses to set the cell division times. It then iterates over microwells and over time and calculates the number of cells and the number of differentiated cells. Full details are given in the paper.

This function simulates cell proliferation and differentiation in microwells over T hours. Each simulated microwell starts with a different number of cells (1 - 10).
The cells can divide or die at rates which were obtained by single cell experiments.
The cells can differentiate either in a constant rate or in a rate that logistically depends on the instantaneous number of cells in a microwell.

Function inputs:
trial - Number of simulates microwells. default is 100
T     - maximal simulated hours. default is 96
DD    - division destiny. a division density can be assigned to the cells which will limit the number of divisions. if it is 0 (default) then cells have unlimited divisions
delayall   - delay the division time of all cells. default is 1 (no delay)
delaydiff  - delay division time of differentiated cells. default is 1 (no delay)
isLogistic - simulate either logistic differentiation (1, default) or 0 (constant differentiation)

Function outputs:
TotalCells_OverTime - total number of cells in each time point for each microwell
DiffCells_OverTime  - total number of differentiated cells in each time point for each microwell

Parameters uploaded from the 'parameter_mat':
logistic differentiation parameters: maximum, Nc, rate
tempbeta - distribution used to add noise to the simulation
weil1d and weil2d - Weibull distributions used to draw first and subsequent division times.
