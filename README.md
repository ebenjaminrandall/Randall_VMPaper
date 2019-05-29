# Randall_VMPaper
Supplemental material 

This repository contains Supplemental material for the manuscript "A model-based analysis of autonomic function in response to the Valsalva maneuver". The SupplementalMaterials.pdf file contains the data (ECG, blood pressure, heart rate, and intrathoracic pressure) for the 34 control subjects (128 individual datasets) and 5 patients (9 individual datasets).

Also included is the associated code for a representative control subject with the optimized parameter values. The driver is the DriverBasic.m file which calls load_global.m (parameter values) and model_sol.m, which calls the Fortran compiled executable driver. The model is located in the file driver_baroreflex.f file, and all other associated Fortran codes are compiled into the driver executable with the command 

gfortran -o driver *.f.

For any questions, please feel free to email ebenjaminrandall@gmail.com. 
