# Predictive_Potential_Temperature_Indices
Scripts for generating results and figures for determining predictive potential of first/second general temperature index for thermodynamic properties.

**File Name	Description**

**Combined_corr_curves.m**	Generates combined correlation curves plots, grouped by experimental data, for the 2 pairs of data, highlighting the curves' peaks and intersection between two curves with highest correlation in some given alpha.

**Good_alpha_intervals.m**	Generates 2 correlation curve plots, one for each pair of data, highlighting the alpha intervals which results in good correlation.

**Scatter_plots.m**	Generates 2 scatter plots, one for each pair of data.

**GoldenSectionSearch_Maximum.m**	IMPORTANT: This must be placed in the working directory as it is referenced in the three other scripts. This is the Octave 7.2 implementation of the Golden Section Search algorithm. I manually translated this from the Python code. 
