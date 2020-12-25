# Galaxy-SED-Fitting-MCMC
This code fits the spectral energy distribution of a given galaxy with a multiple bursts composite stellar population models (Bruzual and Charlot 2003) using Monte Carlo Markov Chain algorithm.
The Swift + MaNGA Catalog consists of a sample of 150 galaxies. The file swim_ba_agn_cut.txt contains a sub-sample of 90 galaxies selected on the basis of axis ratio (>0.6) and non-AGNs.
The mcmc algorithm lies in main_mcmc.py which runs on the sub sample of 90 galaxies. Line 40 and 41 of main_mcmc.py specifies the index of the galaxy on which the mcmc is performed. The file 
Data.zip consists of two galaxies with MaNGA ID 1-23890 (extremely quenched) and 1-37167 (highly star forming). To run the mcmc on 1-23890 put i=0,1 at line 40 and 41 of main_mcmc.py and for 1-37167
put i=4,5 at line 40 and 41 of main_mcmc.py.
