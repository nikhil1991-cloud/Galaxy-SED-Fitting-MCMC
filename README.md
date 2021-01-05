# Galaxy-SED-Fitting-MCMC
This code fits the spectral energy distribution of a given galaxy with a multiple bursts composite stellar population models (Bruzual and Charlot 2003) using Monte Carlo Markov Chain algorithm.
The Swift + MaNGA Catalog consists of a sample of 150 galaxies. The file swim_ba_agn_cut.txt contains a sub-sample of 90 galaxies selected on the basis of axis ratio (>0.6) and non-AGNs.
The mcmc algorithm lies in main_mcmc.py which runs on the sub sample of 90 galaxies. The SFH is derived using code sfr.py which is at line 117 (for initial values of E(B-V) and [Fe/H]) and line 145 (and for random walks in E(B-V) vs [Fe/H] space) of main_mcmc.py. The code sfr.py uses functions from the codes Matrix.py, generate_sigma.py and generate_alphas.py. Once the SFHs and model observables are derived for a specific value of E(B-V) and [Fe/H], we use the function Chi_stats.py to calculate the chi-square between data and observables at lines 122 and 149 of main_mcmc.py.

Line 40 and 41 of main_mcmc.py specifies the index of the galaxy on which the mcmc is performed. The file 
Data.zip consists of two galaxies with MaNGA ID 1-23890 (extremely quenched) and 1-37167 (highly star forming). They are named as SwiM_(MaNGA ID).fits. The information about the file can be found in the header. 
To run the mcmc on 1-23890 put i=0,1 at line 40 and 41 of main_mcmc.py and for 1-37167
put i=4,5 at line 40 and 41 of main_mcmc.py. Similarly, for each of these galaxies,depending on their redshifts, we have generated the bc03 model grids for 30 Ages, which are stored in the file bc03_model_grid_z.zip. They are named as bc03_ssp_(the redshift of that galaxy).fits. For a specific galaxy specify the bc03 model grid depending on its redshift at line 48.
