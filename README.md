# Randall_VMPaper
Supplemental material 

This repository contains Supplemental material for the manuscript "A model-based analysis of autonomic function in response to the Valsalva maneuver". The pdf files entitled SupplementalMaterials_ contain the figures for the data (ECG, blood pressure, heart rate, and intrathoracic pressure) and model predictions for the 35 control subjects (128 individual datasets) and 5 patients (7 individual datasets). SupplementalMaterials_Tables.pdf contains the tables of the descriptions for the control subjects, patient data and clinical ratios for each individual control subject dataset and patient dataset, optimized parameter values for each individual control subject dataset and patient dataset, and results from test for identifiability of parameter subset. 

Also included is the associated code for a representative control subject with the optimized parameter values. The driver is the DriverBasic.m file which calls load_global.m (parameter values) and model_sol.m (solve model), which calls the Fortran compiled executable driver. The model is in the driver_baroreflex.f file, and all other associated Fortran codes are compiled into the driver executable with the command in the command line 

gfortran -o driver *.f.

For any questions, please feel free to email ebenjaminrandall@gmail.com. 
