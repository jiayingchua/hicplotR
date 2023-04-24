# hicplotR

This repository was created for Cranfield BIX Group Project, Group 2.
It contains the R script that performs HiC analysis on the server. The code is run via command line.
The script that is currently on the server, and that the application runs, is PlottingHiC_v3.R

The code utilises the strawr (https://github.com/aidenlab/straw/tree/master/R) and plotgardener (https://phanstiellab.github.io/plotgardener/) to read and plot the Hi-C contact maps from the .hic file.
A separately created function reads the .assembly file to plot the black boxes along the contact map.
Outputs are exported as an image (.png) to a specifed folder.

Jia Ying Chua
April 2023
