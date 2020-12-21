import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import math
import pandas as pd
import time
import random
from scipy.linalg import solve
from scipy import optimize
from scipy.linalg import block_diag
import seaborn as sns
import os
os.chdir('/Users/nikhil/code/Swiftmatch/MPL-7/MCMC_BC03/Chi-sq-mcmc')
from Matrix import Lick_matrix,D4000_matrix,Color_matrix
from Chi_stats import chi_square,chi_square_color
from Generate_sigma import gen_lick_sigma,gen_D4000_sigma,gen_color_sigma
from Generate_alphas import generate_alphas

#Define some constants
fAB = 3.631*math.pow(10,-20)
c = 3*math.pow(10,10)
Tprime = math.pow(10,9)
alpha_min = 1e-5
alpha_max = None
minE=0
maxE=1.5
minZ=-2.23
maxZ=0.5
time_start = time.clock()

#Read Lick indices names
mf = pd.read_csv('/Users/nikhil/code/Newtext/Matchtxt/spec1.txt', comment='#', header=None, delim_whitespace=True)
Indices = np.array(mf[2])

#specify path to text file that contains MaNGA IDs of galaxies with AGN cut and axis ratio cut
with open('/Users/nikhil/code/Newtext/Matchtxt/swim_ba_agn_cut.txt') as f:
     Line0 = [line.rstrip('\n') for line in open('/Users/nikhil/code/Newtext/Matchtxt/swim_ba_agn_cut.txt')]
     
#Read each galaxy data
i=0
for i in range (0,1):#np.shape(Line0)[0]):
 #get galaxy info from drpall
 drpall = fits.open('/Users/nikhil/Data/MaNGAPipe3D/Newmanga/drpall-v2_3_1.fits')
 tbdata = drpall[1].data
 indx = np.where(tbdata['mangaid'] == Line0[i])
 zeta = tbdata['nsa_z'][indx][0] #Redshift
 #Open the pre computed SSP integrated observables for the specific redshift.
 DATA = fits.open('/Users/nikhil/Model_bc03/ALLMETAL/bc03_ssp_'+str(zeta)+'_.fits')
 #Select ages index you wish to include
 t_index = np.array([1,5,8,11,14,17,20,24,27,29])
 Ages = DATA['Ages'].data[t_index]
 #Specify Lick Red, Lick Blue, Lick I, D4000_Red, D4000_Blue, E(B-V), [Fe/H] from DATA
 LR = DATA[8].data
 LB = DATA[9].data
 LI = DATA[10].data
 DI = DATA[11].data
 EBV = DATA[13].data
 Metal = DATA[12].data
 Delta_metal = Metal[3] - Metal[2]
 Delta_Ebv = EBV[3] - EBV[2]
 SD4 = DATA[3].data[:,:,:,t_index]
 SPR = DATA[1].data[:,:,:,t_index]
 SPB = DATA[2].data[:,:,:,t_index]
 SPI = DATA[0].data[:,:,:,t_index]
 SMg = DATA[4].data[:,:,:,t_index]
 K_lick = (DATA['lmid_I'].data - DATA['lmid_B'].data)/(DATA['lmid_R'].data - DATA['lmid_B'].data)
 Delta_lmd = DATA['Dlambda_I'].data
 SMgAB = DATA[5].data#AB Mag normalization
 
 
 
 #Calculate the time intervals for the SFH
 T_age = np.log10(Ages)
 Tmid = (T_age[1:] + T_age[:-1]) / 2
 Tbeg = np.zeros(np.shape(T_age))
 Tend = np.zeros(np.shape(T_age))
 Tbeg[0:len(Ages)-1] = Tmid
 Tbeg[-1] = np.log10(14)
 Tend[1:len(Ages)] = Tmid
 Tend[0] = -9
 delta_T = 10**Tbeg - 10**Tend
 #Open the data file
 #hdu = fits.open("/Volumes/Nikhil/MPL-7_Files/SwiM_binned/SwiM_"+str(Line0[i])+".fits")
 hdu = fits.open("/Users/nikhil/Data/Test/SwiM_"+str(Line0[i])+".fits")
 #We are dealing with integrated bin within 1.5 Re so the integrated binned data is stored as the last element
 #Read mags from data
 Nanomags_data = hdu[22].data[:,-1]
 Mags_data = 22.5-2.5*np.log10(Nanomags_data)
 Nanomags_error = np.sqrt(hdu[23].data[:,-1])
 Mags_error = (2.5)*(Nanomags_error/Nanomags_data)
 #Read Lick indices and d4000 from data
 Dn4000_data = hdu[18].data[43][-1]
 Dn4000_error = hdu[19].data[43][-1]
 Lick_indices_data = hdu[18].data[0:42,:][:,-1]
 Lick_indices_error = hdu[19].data[0:42][:,-1]
 
 #generate alpha matrices
 Lick_list = [3,4,13,14,21,22,23,24] #fe4383,fe5270,fe5335,HdA
 m1_list,m2_list = [0,1,4,4,3,3],[3,3,6,5,6,7]  #w2-u,w1-u,g-i,g-r,u-i,u-z
 mags_info = np.array(['W2','W1','M2','u','g','r','i','z'])
 
 P=2000#Specify epochs
 E_step = np.array([0.05])#set step size for E(B-V)
 Z_step = np.array([0.05])#set step size for [Fe/H]
 
 #Set initial E and Z
 Initial_E = 1.5
 Initial_Z = 0.2
 
 #dof for regression and mcmc
 dof_mcmc = len(m1_list) + len(Lick_list)
 dof_lr = len(m1_list) + len(Lick_list) + 1 - len(Ages)
 
 #compute model observables for given Z and E using code sfr.py
 Given_Z = Initial_Z
 Given_E = Initial_E
 exec(open("/Users/nikhil/code/Swiftmatch/MPL-7/MCMC_BC03/Chi-sq-mcmc/sfr.py").read())
 Model_D4000_I = D4000_Model
 Model_LICK_I = LICK_Model
 Model_MAGS_I = MAGS_Model
 #compute initial chi-square
 CHI_initial = np.sum((chi_square(Dn4000_data,Model_D4000_I,Dn4000_error),np.sum(chi_square(Lick_indices_data,Model_LICK_I,Lick_indices_error)[Lick_list]),np.sum(chi_square_color(m1_list,m2_list,Mags_data,Model_MAGS_I,Mags_error))))
 #Initialize all arrays with the initial observables and parameters
 output_parameters = [[CHI_initial,Initial_Z,Initial_E,Model_D4000_I,Model_LICK_I[4],Model_LICK_I[13],Model_LICK_I[14],Model_LICK_I[21],Model_MAGS_I[0]-Model_MAGS_I[3],Model_MAGS_I[1]-Model_MAGS_I[3],Model_MAGS_I[4]-Model_MAGS_I[6],Model_MAGS_I[4]-Model_MAGS_I[5],Model_MAGS_I[3]-Model_MAGS_I[6],Model_MAGS_I[3]-Model_MAGS_I[7],alpha_bestfit[0],alpha_bestfit[1],alpha_bestfit[2],alpha_bestfit[3],alpha_bestfit[4],alpha_bestfit[5],alpha_bestfit[6],alpha_bestfit[7],alpha_bestfit[8],alpha_bestfit[9]]]
 CHI_Prop = [CHI_initial]
 #Generate array that lists the names of all the output parameters called 'output_info'
 output_info = ['']*(np.shape(output_parameters)[1])
 output_info[0:4] = ['Chisq','[Fe/H]','E(B-V)','Dn4000']
 output_info[4:8] = Indices[[4,13,14,21]]
 output_info[8:14] = np.char.add(mags_info[m1_list],mags_info[m2_list])
 output_info[14:] = ['a0','a1','a2','a3','a4','a5','a6','a7','a8','a9']

 #Start walking
 m=0
 for m in range (0,P):
           #set step size
           dE = E_step[-1]*np.random.normal(scale=1.0)
           dZ = Z_step[-1]*np.random.normal(scale=1.0)
           #New parameters
           E_trial = np.clip(output_parameters[-1][2] + dE,minE,maxE)
           Z_trial = np.clip(output_parameters[-1][1] + dZ,minZ,maxZ)
           #Compute chi-sq for new parameters using 'step.py'
           Given_E = E_trial
           Given_Z = Z_trial
           exec(open("/Users/nikhil/code/Swiftmatch/MPL-7/MCMC_BC03/Chi-sq-mcmc/sfr.py").read())
           Model_D4000 = D4000_Model
           Model_LICK = LICK_Model
           Model_MAGS = MAGS_Model
           CHI = np.sum((chi_square(Dn4000_data,Model_D4000,Dn4000_error),np.sum(chi_square(Lick_indices_data,Model_LICK,Lick_indices_error)[Lick_list]),np.sum(chi_square_color(m1_list,m2_list,Mags_data,Model_MAGS,Mags_error))))
           CHI_Prop.append(CHI)
           #Accept
           if CHI<output_parameters[-1][0]:
              output_parameters.append([CHI,Z_trial,E_trial,Model_D4000,Model_LICK[4],Model_LICK[13],Model_LICK[14],Model_LICK[21],Model_MAGS[0]-Model_MAGS[3],Model_MAGS[1]-Model_MAGS[3],Model_MAGS[4]-Model_MAGS[6],Model_MAGS[4]-Model_MAGS[5],Model_MAGS[3]-Model_MAGS[6],Model_MAGS[3]-Model_MAGS[7],alpha_bestfit[0],alpha_bestfit[1],alpha_bestfit[2],alpha_bestfit[3],alpha_bestfit[4],alpha_bestfit[5],alpha_bestfit[6],alpha_bestfit[7],alpha_bestfit[8],alpha_bestfit[9]])
           else:
              prob = np.exp(-(CHI-output_parameters[-1][0]))
              rand = np.random.rand(1)
              if rand <= prob:
                 output_parameters.append([CHI,Z_trial,E_trial,Model_D4000,Model_LICK[4],Model_LICK[13],Model_LICK[14],Model_LICK[21],Model_MAGS[0]-Model_MAGS[3],Model_MAGS[1]-Model_MAGS[3],Model_MAGS[4]-Model_MAGS[6],Model_MAGS[4]-Model_MAGS[5],Model_MAGS[3]-Model_MAGS[6],Model_MAGS[3]-Model_MAGS[7],alpha_bestfit[0],alpha_bestfit[1],alpha_bestfit[2],alpha_bestfit[3],alpha_bestfit[4],alpha_bestfit[5],alpha_bestfit[6],alpha_bestfit[7],alpha_bestfit[8],alpha_bestfit[9]])
           #acceptance rate check and adjustment of step size
           if m%100==0:
           
              AR = (np.shape(output_parameters)[0]/np.shape(CHI_Prop)[0])*100


              if AR > 30:
                 E_step = np.append(E_step,E_step[-1]*2)
                 Z_step = np.append(Z_step,Z_step[-1]*2)

              if AR < 15:
                 E_step = np.append(E_step,E_step[-1]/2)
                 Z_step = np.append(Z_step,Z_step[-1]/2)


              
              
 print (time.clock() - time_start, "seconds")
 #Store all values
 OUTPUT = np.array(output_parameters)
 OUTPUT[:,0] = OUTPUT[:,0]/dof_mcmc #Reduced chi-square

 plt.subplot(2,2,1)
 plt.scatter(OUTPUT[:,3],OUTPUT[:,7],c=OUTPUT[:,0],s=20,vmin=OUTPUT[:,0].min(),vmax=OUTPUT[:,0].max())
 plt.scatter(Dn4000_data,Lick_indices_data[21],c='k',s=35)
 plt.xlabel('Dn4000')
 plt.ylabel('HdA')
 
 plt.subplot(2,2,2)
 plt.scatter(OUTPUT[:,8],OUTPUT[:,10],c=OUTPUT[:,0],s=20,vmin=OUTPUT[:,0].min(),vmax=OUTPUT[:,0].max())
 plt.scatter(Model_MAGS[0]-Model_MAGS[3],Model_MAGS[4]-Model_MAGS[6],c='k',s=35)
 plt.xlabel('w2-u')
 plt.ylabel('g-i')
 
 plt.subplot(2,2,3)
 plt.scatter(OUTPUT[:,3],OUTPUT[:,4],c=OUTPUT[:,0],s=20,vmin=OUTPUT[:,0].min(),vmax=OUTPUT[:,0].max())
 plt.scatter(Dn4000_data,Lick_indices_data[3],c='k',s=35)
 plt.xlabel('Dn4000')
 plt.ylabel('Fe4383')
 
 plt.subplot(2,2,4)
 plt.scatter(OUTPUT[:,1],OUTPUT[:,2],c=OUTPUT[:,0],s=20,vmin=OUTPUT[:,0].min(),vmax=OUTPUT[:,0].max())
 plt.xlabel('[Fe/H]')
 plt.ylabel('E(B-V)')

 #Generate a header file
 #hdu0 = fits.PrimaryHDU(OUTPUT)
 #new_hdul = fits.HDUList([hdu0])
 #new_hdul.writeto("/Volumes/Nikhil/MCMC_New/New_sampling/MCMC_Integrated/"+str(Line0[i])+"-STATS_1bin.fits")
