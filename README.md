# AuNP_tilt_cryoem
Code for reconstruction of Gold nanoparticle 3D coordinates from Cryo-EM tilt-pairs

This is an implementation of a method currently under review for publication.
It has been tested in MATLAB R2020b and requires the image processing toolbox.
We have provided the real and simulated cryo-EM datasets on which the code can be run.

The programs should be run in the following sequential order to recover 3D coordinates for the provided real tilt-pairs:
AuNP_reconstruction.m -> initial_aggregation.m -> final_aggregation.m -> pick_correspond.m -> final_aggregation.m

To choose which dataset to use, set the variable 'real' in AuNP_reconstruction.m: 0 for simulated data, 1 for real data.
When using simulated data, the run of initial_aggregation.m can be skipped since the reference model coordinates are
used for alignment in final_aggregation.m in this case.

Additionally, circle_plots.m plots the final list of 3D structures while save_struct.m saves the data.
find_Au.m, plot_mols.m, and kabsch.m are helper functions.
make_real_pdbs.py and golds.pdb are used to convert saved structures into pdb files with all gold atoms included.

The functionality of each program is described with in-line comments.
