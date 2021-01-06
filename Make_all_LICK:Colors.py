import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import scipy.stats as st
import math
from scipy import signal
from astropy import wcs
import sys
from scipy import interpolate
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import scipy.stats as st
from scipy import signal
from astropy import wcs
import sys
from scipy import interpolate
import matplotlib.patches as mpatches
from skimage.transform import rotate
from scipy.stats import norm
import scipy.stats as stats
from numpy import inf
from astropy.cosmology import WMAP9 as cosmo
from astropy.stats import sigma_clip
from astropy.stats import biweight_location
import pandas as pd
from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord
import extinction
from random import randint
import itertools
import ezgal
from extinction import apply
import time

time_start = time.clock()

hdu = fits.open("/Users/nikhil/Data/MaNGAPipe3D/manga-7495-12703-LOGCUBE.fits")
Waveh = hdu['WAVE'].data

with open('/Users/nikhil/code/Newtext/Matchtxt/SwiM_ID.txt') as f:
    Line0 = [line.rstrip('\n') for line in open('/Users/nikhil/code/Newtext/Matchtxt/Swim_ID.txt')]
IZ =np.zeros(np.shape(Line0))
q=0
for q in range (0,np.shape(Line0)[0]):
  ID = Line0[q]
  drpall = fits.open('/Users/nikhil/Data/MaNGAPipe3D/Newmanga/drpall-v2_3_1.fits')
  tbdata = drpall[1].data
  indx = np.where(tbdata['mangaid'] == Line0[q])
  IZ[q] = tbdata['nsa_z'][indx][0]
  
L_solar = 3.839*math.pow(10,33)
fAB = 3.631*math.pow(10,-20)
c = 2.99792*math.pow(10,10)
#Read photometric filters
mFUV = pd.read_csv('/Users/nikhil/code/Newtext/Filt/FUV_Prime.txt', comment='#', header=None, delim_whitespace=True)
mNUV = pd.read_csv('/Users/nikhil/code/Newtext/Filt/NUV_Prime.txt', comment='#', header=None, delim_whitespace=True)
mW1  = pd.read_csv('/Users/nikhil/code/Newtext/Filt/W1_Prime.txt', comment='#', header=None, delim_whitespace=True)
mW2  = pd.read_csv('/Users/nikhil/code/Newtext/Filt/W2_Prime.txt', comment='#', header=None, delim_whitespace=True)
mM2  = pd.read_csv('/Users/nikhil/code/Newtext/Filt/M2_Prime.txt', comment='#', header=None, delim_whitespace=True)
mU  = pd.read_csv('/Users/nikhil/code/Newtext/Filt/U_Prime.txt', comment='#', header=None, delim_whitespace=True)
mG  = pd.read_csv('/Users/nikhil/code/Newtext/Filt/G_Prime.txt', comment='#', header=None, delim_whitespace=True)
mR  = pd.read_csv('/Users/nikhil/code/Newtext/Filt/R_Prime.txt', comment='#', header=None, delim_whitespace=True)
mI  = pd.read_csv('/Users/nikhil/code/Newtext/Filt/I_Prime.txt', comment='#', header=None, delim_whitespace=True)
mZ  = pd.read_csv('/Users/nikhil/code/Newtext/Filt/Z_Prime.txt', comment='#', header=None, delim_whitespace=True)
#Read Data of photometric filters
RFUV = np.array(mFUV)
RNUV = np.array(mNUV)
RW1 = np.array(mW1)
RW2 = np.array(mW2)
RM2 = np.array(mM2)
RU = np.array(mU)
RG = np.array(mG)
RR = np.array(mR)
RI = np.array(mI)
RZ = np.array(mZ)
#Read Fluxes
RfFUV = RFUV[:,1]
RfNUV = RNUV[:,1]
RfW1 = RW1[:,1]
RfW2 = RW2[:,1]
RfM2 = RM2[:,1]
RfU = RU[:,1]
RfG = RG[:,1]
RfR = RR[:,1]
RfI = RI[:,1]
RfZ = RZ[:,1]
#Read Wavelengths
RwFUV =RFUV[:,0]
RwNUV = RNUV[:,0]
RwW1 = RW1[:,0]
RwW2 = RW2[:,0]
RwM2 = RM2[:,0]
RwU = RU[:,0]
RwG = RG[:,0]
RwR = RR[:,0]
RwI = RI[:,0]
RwZ = RZ[:,0]
#Read Lick indices
mf = pd.read_csv('/Users/nikhil/code/Newtext/Matchtxt/spec1.txt', comment='#', header=None, delim_whitespace=True)
Indices = np.array(mf[2])
H_Blue = np.array(mf[9])
L_Blue = np.array(mf[8])
H_Red = np.array(mf[13])
L_Red = np.array(mf[12])
H_Index = np.array(mf[5])
L_Index = np.array(mf[4])
Unit_indx = np.array(mf[16])
KJ = np.where(Unit_indx == 'mag')
KL = Unit_indx*0
KL[KJ] = 1
Rmid = (H_Red + L_Red)/2
Bmid = (H_Blue + L_Blue)/2
Imid = (H_Index + L_Index)/2
Delt_In = H_Index - L_Index
Delt_R = H_Red - L_Red
Delt_B = H_Blue - L_Blue
#Initiate Arrays for lick
Ages = np.array([0.0030,0.0035,0.0054,0.0063,0.0099,0.0114,0.0181,0.0206,0.0331,0.0374,0.0605,0.0676,0.1103,0.1224,0.2012,0.2213,0.3670,0.4003,0.6692,0.7238,1.2205,1.3090,2.2258,2.3673,4.0592,4.2810,7.4026,7.7417,13.500,14.000])
IK_prime = Ages
EBV = np.linspace(0,1.5,30)
C_red = np.zeros((np.shape(mf)[0],np.shape(EBV)[0],np.shape(IK_prime)[0]))
C_blue = np.zeros((np.shape(mf)[0],np.shape(EBV)[0],np.shape(IK_prime)[0]))
C_index = np.zeros((np.shape(mf)[0],np.shape(EBV)[0],np.shape(IK_prime)[0]))
Continuum = np.zeros((np.shape(mf)[0],np.shape(EBV)[0],np.shape(IK_prime)[0]))
LICK = np.zeros((np.shape(mf)[0],np.shape(EBV)[0],np.shape(IK_prime)[0]))
Photo_index = np.zeros((10,np.shape(EBV)[0],np.shape(IK_prime)[0]))
PhotoAB_index = np.zeros((10,np.shape(EBV)[0]))
llmid_R = np.zeros((np.shape(mf))[0])
llmid_B = np.zeros((np.shape(mf))[0])
llmid_I = np.zeros((np.shape(mf))[0])
Dlambda_I = np.zeros((np.shape(mf))[0])
D4= np.zeros((2,np.shape(EBV)[0],np.shape(IK_prime)[0]))
D4_nu= np.zeros((2,np.shape(EBV)[0],np.shape(IK_prime)[0]))
W = np.array([-2.2490,-1.6464,-0.6392,-0.3300,0.0932,0.5595])
red_shift = np.sort(IZ)
e=0
for e in range (0,np.shape(W)[0]):
 z=0
 for z in range (0,np.shape(red_shift)[0]):
  s=0
  for s in range (0,np.shape(EBV)[0]):
    Ext = EBV[s]
    Av = 3.1 * Ext
    Frame = pd.read_csv('/Users/nikhil/bc03_ISEDs/bc2003_hr_'+str(W[e])+'.44', comment='#', header=None, delim_whitespace=True)
    FrameA = np.array(Frame)
    ll = FrameA[:,0]
    Data = FrameA[:,1:]*L_solar
    ll_shift = ll*(1+red_shift[z])
    Alambda = extinction.calzetti00(ll, Av, 3.1)
    Alambda_red = extinction.calzetti00(ll_shift, Av, 3.1)
    Data_fl = np.zeros(np.shape(Data))
    Data_O = np.zeros(np.shape(Data))
    w=0
    for w in range (0,np.shape(Data)[1]):
        Data_O[:,w] = extinction.apply(Alambda_red,Data[:,w])
        Data_fl[:,w] = extinction.apply(Alambda,Data[:,w])
    #Color computation in observed frame
    lmd_shift =ll_shift
    FL_red = Data_O/(1+red_shift[z])
    #Color_computation for ssp ages
    #Read fluxes for broad bands
    frFUV_red = FL_red[np.where((lmd_shift>=RwFUV.min())&(lmd_shift<=RwFUV.max())),:][0]
    frNUV_red = FL_red[np.where((lmd_shift>=RwNUV.min())&(lmd_shift<=RwNUV.max())),:][0]
    frW1_red = FL_red[np.where((lmd_shift>=RwW1.min())&(lmd_shift<=RwW1.max())),:][0]
    frW2_red = FL_red[np.where((lmd_shift>=RwW2.min())&(lmd_shift<=RwW2.max())),:][0]
    frM2_red = FL_red[np.where((lmd_shift>=RwM2.min())&(lmd_shift<=RwM2.max())),:][0]
    frU_red = FL_red[np.where((lmd_shift>=RwU.min())&(lmd_shift<=RwU.max())),:][0]
    frG_red = FL_red[np.where((lmd_shift>=RwG.min())&(lmd_shift<=RwG.max())),:][0]
    frR_red = FL_red[np.where((lmd_shift>=RwR.min())&(lmd_shift<=RwR.max())),:][0]
    frI_red = FL_red[np.where((lmd_shift>=RwI.min())&(lmd_shift<=RwI.max())),:][0]
    frZ_red = FL_red[np.where((lmd_shift>=RwZ.min())&(lmd_shift<=RwZ.max())),:][0]
    #Response wave lengths for broad bands
    RwFUV_red = lmd_shift[np.where((lmd_shift>=RwFUV.min())&(lmd_shift<=RwFUV.max()))]
    RwNUV_red = lmd_shift[np.where((lmd_shift>=RwNUV.min())&(lmd_shift<=RwNUV.max()))]
    RwW1_red = lmd_shift[np.where((lmd_shift>=RwW1.min())&(lmd_shift<=RwW1.max()))]
    RwW2_red = lmd_shift[np.where((lmd_shift>=RwW2.min())&(lmd_shift<=RwW2.max()))]
    RwM2_red = lmd_shift[np.where((lmd_shift>=RwM2.min())&(lmd_shift<=RwM2.max()))]
    RwU_red = lmd_shift[np.where((lmd_shift>=RwU.min())&(lmd_shift<=RwU.max()))]
    RwG_red = lmd_shift[np.where((lmd_shift>=RwG.min())&(lmd_shift<=RwG.max()))]
    RwR_red = lmd_shift[np.where((lmd_shift>=RwR.min())&(lmd_shift<=RwR.max()))]
    RwI_red = lmd_shift[np.where((lmd_shift>=RwI.min())&(lmd_shift<=RwI.max()))]
    RwZ_red = lmd_shift[np.where((lmd_shift>=RwZ.min())&(lmd_shift<=RwZ.max()))]
    #Interpolating responses at response wavelengths of broadbands
    f_FUV = interpolate.interp1d(RwFUV,RfFUV)
    f_NUV = interpolate.interp1d(RwNUV,RfNUV)
    f_W1 = interpolate.interp1d(RwW1,RfW1)
    f_W2 = interpolate.interp1d(RwW2,RfW2)
    f_M2 = interpolate.interp1d(RwM2,RfM2)
    f_U = interpolate.interp1d(RwU,RfU)
    f_G = interpolate.interp1d(RwG,RfG)
    f_R = interpolate.interp1d(RwR,RfR)
    f_I = interpolate.interp1d(RwI,RfI)
    f_Z = interpolate.interp1d(RwZ,RfZ)
    #Calculate the response at response wavelengths
    RfFUV_red = f_FUV(RwFUV_red)
    RfNUV_red = f_NUV(RwNUV_red)
    RfW1_red = f_W1(RwW1_red)
    RfW2_red = f_W2(RwW2_red)
    RfM2_red = f_M2(RwM2_red)
    RfU_red = f_U(RwU_red)
    RfG_red = f_G(RwG_red)
    RfR_red = f_R(RwR_red)
    RfI_red = f_I(RwI_red)
    RfZ_red = f_Z(RwZ_red)
    
    #plt.plot(lmd_shift,np.log10(FL_red[:,16]),c='r')
    #plt.plot(ll,np.log10(Data_O[:,16]),c='y')
    #plt.plot(RwG_red,np.log10(frG_red[:,16]),c='g')
    #plt.plot(RwG,np.log10(math.pow(10,32)*RfG),c='lightgreen')
    #plt.xlim([3500,6000])
    #plt.ylim([30,32])

    #Convolve with response (lambda*f_lambda*R_P)
    PFUV = np.transpose(frFUV_red)*RfFUV_red*RwFUV_red
    PNUV = np.transpose(frNUV_red)*RfNUV_red*RwNUV_red
    PW1 = np.transpose(frW1_red)*RfW1_red*RwW1_red
    PW2 = np.transpose(frW2_red)*RfW2_red*RwW2_red
    PM2 = np.transpose(frM2_red)*RfM2_red*RwM2_red
    PU = np.transpose(frU_red)*RfU_red*RwU_red
    PG = np.transpose(frG_red)*RfG_red*RwG_red
    PR = np.transpose(frR_red)*RfR_red*RwR_red
    PI = np.transpose(frI_red)*RfI_red*RwI_red
    PZ = np.transpose(frZ_red)*RfZ_red*RwZ_red
    #Fluxes for AB in Flambda (lambda*f_AB_lambda*R_P)
    PFUV1 = (fAB*c*(RfFUV_red/(np.power(RwFUV_red,2))))*RwFUV_red
    PNUV1 = (fAB*c*(RfNUV_red/(np.power(RwNUV_red,2))))*RwNUV_red
    PW11 = (fAB*c*(RfW1_red/(np.power(RwW1_red,2))))*RwW1_red
    PW21 = (fAB*c*(RfW2_red/(np.power(RwW2_red,2))))*RwW2_red
    PM21 = (fAB*c*(RfM2_red/(np.power(RwM2_red,2))))*RwM2_red
    PU1 = (fAB*c*(RfU_red/(np.power(RwU_red,2))))*RwU_red
    PG1 = (fAB*c*(RfG_red/(np.power(RwG_red,2))))*RwG_red
    PR1 = (fAB*c*(RfR_red/(np.power(RwR_red,2))))*RwR_red
    PI1 = (fAB*c*(RfI_red/(np.power(RwI_red,2))))*RwI_red
    PZ1 = (fAB*c*(RfZ_red/(np.power(RwZ_red,2))))*RwZ_red
    t0=0
    for t0 in range (0,np.shape(IK_prime)[0]):
        #Numerator of equation 7 in Blanton et al.
        Photo_index[0,s,t0] = np.trapz(PFUV[t0,:],x=RwFUV_red)
        Photo_index[1,s,t0] = np.trapz(PNUV[t0,:],x=RwNUV_red)
        Photo_index[2,s,t0] = np.trapz(PW1[t0,:],x=RwW1_red)
        Photo_index[3,s,t0] = np.trapz(PW2[t0,:],x=RwW2_red)
        Photo_index[4,s,t0] = np.trapz(PM2[t0,:],x=RwM2_red)
        Photo_index[5,s,t0] = np.trapz(PU[t0,:],x=RwU_red)
        Photo_index[6,s,t0] = np.trapz(PG[t0,:],x=RwG_red)
        Photo_index[7,s,t0] = np.trapz(PR[t0,:],x=RwR_red)
        Photo_index[8,s,t0] = np.trapz(PI[t0,:],x=RwI_red)
        Photo_index[9,s,t0] = np.trapz(PZ[t0,:],x=RwZ_red)
        #Denominator of equation 7 in Blanton et al.
        PhotoAB_index[0,s] = np.trapz(PFUV1,x=RwFUV_red)
        PhotoAB_index[1,s] = np.trapz(PNUV1,x=RwNUV_red)
        PhotoAB_index[2,s] = np.trapz(PW11,x=RwW1_red)
        PhotoAB_index[3,s] = np.trapz(PW21,x=RwW2_red)
        PhotoAB_index[4,s] = np.trapz(PM21,x=RwM2_red)
        PhotoAB_index[5,s] = np.trapz(PU1,x=RwU_red)
        PhotoAB_index[6,s] = np.trapz(PG1,x=RwG_red)
        PhotoAB_index[7,s] = np.trapz(PR1,x=RwR_red)
        PhotoAB_index[8,s] = np.trapz(PI1,x=RwI_red)
        PhotoAB_index[9,s] = np.trapz(PZ1,x=RwZ_red)
        
        
        finterp_new = interpolate.interp1d(ll,Data_fl[:,t0],kind='linear')
        New_F = finterp_new(Waveh)#convert to fnu
        #Dn4000_Computation
        H1 = np.where(Waveh < 4000)
        L1 = np.where(Waveh > 4100)
        H10 = np.max(H1)
        L10 = np.min(L1)
        ar1 = Waveh[H10:L10]
        D4[0,s,t0] = (1/c)*np.trapz(New_F[H10:L10]*(ar1**2),x=ar1)
        H2 = np.where(Waveh < 3850)
        L2 = np.where(Waveh > 3950)
        H20 = np.max(H2)
        L20 = np.min(L2)
        ar2 = Waveh[H20:L20]
        D4_red_blue = np.array((ar1[np.int(np.shape(ar1)[0]/2)],ar2[np.int(np.shape(ar2)[0]/2)]))
        D4[1,s,t0] = (1/c)*np.trapz(New_F[H20:L20]*(ar2**2),x=ar2)
        D4000_C = D4[0,s,:]/D4[1,s,:]
        #LICK computation
        l=0
        for l in range (0,np.shape(mf)[0]):
            High_Blue = np.where(Waveh < L_Blue[l]) #Blue Window
            Low_Blue = np.where(Waveh > H_Blue[l])
            HB0 = np.max(High_Blue)
            LB0 = np.min(Low_Blue)
            DBlue = HB0 - LB0
            Wave_Blue = Waveh[HB0:LB0]
            mid_Blue = np.int(np.shape(Wave_Blue)[0]/2)
            Flux_Blue = New_F[HB0:LB0]
            C_blue[l,s,t0] = np.trapz(Flux_Blue[:],x=Wave_Blue,axis=0)/Delt_B[l]
            High_Red = np.where(Waveh < L_Red[l]) #Red Window
            Low_Red = np.where(Waveh > H_Red[l])
            HR0 = np.max(High_Red)
            LR0 = np.min(Low_Red)
            DRed = HR0 - LR0
            Wave_Red = Waveh[HR0:LR0]
            mid_Red = np.int(np.shape(Wave_Red)[0]/2)
            Flux_Red = New_F[HR0:LR0]
            C_red[l,s,t0] = np.trapz(Flux_Red[:],x=Wave_Red,axis=0)/Delt_R[l]
            High_Index = np.where(Waveh < L_Index[l]) #Index window
            Low_Index = np.where(Waveh > H_Index[l])
            HI0 = np.max(High_Index)
            LI0 = np.min(Low_Index)
            Wave_Index = Waveh[HI0:LI0]
            Dlambda = Wave_Index.max() - Wave_Index.min()
            mid_Index = np.int(np.shape(Wave_Index)[0]/2)
            IndexFlux = New_F[HI0:LI0]
            C_index[l,s,t0] = np.trapz(IndexFlux[:],x=Wave_Index,axis=0)
            llmid_R[l] = Rmid[l]
            llmid_B[l] = Bmid[l]
            llmid_I[l] = Imid[l]
            Dlambda_I[l] = Delt_In[l]
            Continuum[l,s,t0] = ((C_red[l,s,t0] - C_blue[l,s,t0])/(llmid_R[l]-llmid_B[l]))*(llmid_I[l]-llmid_B[l]) + C_blue[l,s,t0]
            LICK[l,s,t0] = Dlambda_I[l] - C_index[l,s,t0]/Continuum[l,s,t0]
    #p=0
    #Maggie = -2.5*np.log10(Photo_index[:,p,-1]/PhotoAB_index[:,p])
    #Mag = Maggie - Maggie[7]
   # D4_T = np.sum(D4[0,p,:])/np.sum(D4[1,p,:])
   # slope = (np.sum(C_red[:,p,:],axis=1) - np.sum(C_blue[:,p,:],axis=1))/(llmid_R-llmid_B)
   # cont = slope*(llmid_I - llmid_B) + np.sum(C_blue[:,p,:],axis=1)
   # Lick = Dlambda_I - (np.sum(C_index[:,p,:],axis=1)/cont)
#time_elapsed = (time.clock() - time_start)
#print('Execution Time='+str(time_elapsed)+'secs')
  hdu0 = fits.PrimaryHDU(C_index)
  hdu1 = fits.ImageHDU(C_red)
  hdu2 = fits.ImageHDU(C_blue)
  hdu3 = fits.ImageHDU(D4)
  hdu4 = fits.ImageHDU(Photo_index)
  hdu5 = fits.ImageHDU(PhotoAB_index)
  hdu6 = fits.ImageHDU(Ages)
  hdu7 = fits.ImageHDU(ll)
  hdu8 = fits.ImageHDU(llmid_R)
  hdu9 = fits.ImageHDU(llmid_B)
  hdu10 = fits.ImageHDU(llmid_I)
  hdu11 = fits.ImageHDU(Dlambda_I)
  new_hdul = fits.HDUList([hdu0,hdu1, hdu2, hdu3, hdu4, hdu5,hdu6,hdu7,hdu8,hdu9,hdu10,hdu11])
  new_hdul.writeto("/Users/nikhil/Model_bc03/ALLMETAL_bc03/bc03_ssp_"+str(red_shift[z])+"_"+str(W[e])+".fits")
    
  hdu = fits.open("/Users/nikhil/Model_bc03/ALLMETAL_bc03/bc03_ssp_"+str(red_shift[z])+"_"+str(W[e])+".fits")

  hdr = hdu[0].header
  hdr['EXTNAME'] = 'C_Index[Index,Ext,Age]'

  hdr1 = hdu[1].header
  hdr1['EXTNAME'] = 'C_Red[Index,Ext,Age]'

  hdr2 = hdu[2].header
  hdr2['EXTNAME'] = 'C_Blue[Index,Ext,Age]'

  hdr3 = hdu[3].header
  hdr3['EXTNAME'] = 'D4000_Red_Blue[Band,Ext,Age]'
  hdr3['Fred_mid'] = D4_red_blue[0]
  hdr3['Fblue_mid'] = D4_red_blue[1]

  hdr4 = hdu[4].header
  hdr4['EXTNAME'] = 'Photo[FNW1W2M2ugriz,Ext,Age]'

  hdr5 = hdu[5].header
  hdr5['EXTNAME'] = 'PhotoAB[FNW1W2M2ugriz,Ext,Age]'
    
  hdr6 = hdu[6].header
  hdr6['EXTNAME'] = 'Ages'

  hdr7 = hdu[7].header
  hdr7['EXTNAME'] = 'ls'
    
  hdr8 = hdu[8].header
  hdr8['EXTNAME'] = 'lmid_R'
    
  hdr9 = hdu[9].header
  hdr9['EXTNAME'] = 'lmid_B'
    
  hdr10 = hdu[10].header
  hdr10['EXTNAME'] = 'lmid_I'
    
  hdr11 = hdu[11].header
  hdr11['EXTNAME'] = 'Dl_I'

    
  hdu.writeto("/Users/nikhil/Model_bc03/ALLMETAL_bc03/bc03_ssp_"+str(red_shift[z])+"_"+str(W[e])+".fits",overwrite=True)

