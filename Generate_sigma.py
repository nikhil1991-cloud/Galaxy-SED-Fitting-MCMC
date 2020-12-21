import numpy as np


"""This code contains functions to generate the elements of Sigma_double_dahsed (equation 14.6). We have m measurements in total which are
   l number of lick indices, 1 Dn4000 and c number of colors such that m = l + 1 + c. For lick indices, the Sigma_double_dahsed matrix elements
   are the first l diagonal matrices of 14.6 = sigma_0,sigma_1,....sigma_l. The next diagonal element is the Sigma_Dn4000 matrix and the last
   diagonal elements are the c Sigma_color matrices. Each of these elements are a square matrix with shape (j,j) where j are the total measured
   parameters (for 10 ages in our case).
"""


def gen_lick_sigma(model_continuum,data_error,lick_index_list,Ages):
    """Generates Sigma_lick for l lick indices defined at j ages.
    
    Input Parameters:
          model_continuum: Array of shape (j,l) (continuum values f_c for SSP ages j)
          
          data_error: Array of shape (1,l)
          
          lick_index_list: Indices of l lick indices selected in the fit amongst all 42 lick indices
          
          Ages: Array of shape (1,j)
    
    Returns:
          
          Sigma_Lick: Array of shape (l,j,j)
    
    """
    model_continuum_shift = np.concatenate((model_continuum[-1,:].reshape(1,-1),model_continuum[0:model_continuum.shape[0]-1,:]),axis=0)
    Continuum = np.repeat(model_continuum_shift[:,np.newaxis,:],np.shape(model_continuum_shift)[0],axis=1)
    Continuum = Continuum[:,:,lick_index_list]
    Lick_error = data_error[[lick_index_list]]
    Sigma_Lick = np.zeros((len(lick_index_list),len(Ages),len(Ages)))
    j=0
    for j in range (0,len(lick_index_list)):
        Sigma_Lick[j,:,:] = (Continuum[:,:,j]*Continuum[:,:,j].T)*(Lick_error[j]**2)
    return Sigma_Lick

def gen_D4000_sigma(Blue_Flux,data_error,Ages):
    """Generate Sigma_Dn4000 for Dn4000 values defined at j ages.
    
    Input parameters:
          Blue_Flux: Array of shape (1,j)
          
          data_error: float64
          
          Ages: Array of shape (1,j)
          
    Returns:
    
          Sigma_Dn4000: Array of shape (1,j,j)
    """
    Blue_Flux_shift = np.roll(Blue_Flux,1)
    Flux_Blue = np.repeat(Blue_Flux_shift[:,np.newaxis],np.shape(Blue_Flux_shift)[0],axis=1)
    Sigma_D4 = (Flux_Blue[:,:]*Flux_Blue[:,:].T)*(data_error**2)
    return Sigma_D4

def gen_color_sigma(l1,l2,Flux_2,data_mag_error,data_magnitudes,Ages):
    """Generates Sigma_color for c colors defined at j ages. The color between broad-band filter 1 and 2 is defined as c=m1-m2 where m1 and m2.
       
    Input Parameters:
    
          l1: List of indices of selected broad-band filter 1 (W1,W2,M2,u,g,r,i,z)
          
          l2: List of indices of selected broad-band filter 2 (W1,W2,M2,u,g,r,i,z)
          
          For example l1 = [0,1,2] and l2 = [3,3,3] will give 3 colors W1-u,W2-u and M2-u
    
          Flux_2: Array of shape (c,j) (Flux in the broad-band filter 2 for SSP ages j)
          
          data_mag_error: Array of shape (1,c) Given by (2.5)*(Nanomags_error/Nanomags_data)
          
          data_magnitudes: Array of shape (1,c) Given by 22.5-2.5*np.log10(Nanomags_data)
          
          Ages: Array of shape (1,j)
    
    Returns:
          
          Sigma_Color: Array of shape (c,j,j)
    
    """
    Flux2 = Flux_2.T
    Flux2_shift = np.concatenate((Flux2[-1,:].reshape(1,-1),Flux2[0:Flux2.shape[0]-1,:]),axis=0)
    F2 = np.repeat(Flux2_shift[:,np.newaxis,:],np.shape(Flux2_shift)[0],axis=1)
    F2 = F2[:,:,l2]
    Color_error = (np.sqrt(data_mag_error[l1]**2 + data_mag_error[l2]**2))
    Data_color = data_magnitudes[l1] - data_magnitudes[l2]
    sigma_cov = 0.84*(10**(-0.84*Data_color))*(Color_error**2)
    Sigma_color = np.zeros((len(l2),len(Ages),len(Ages)))
    j=0
    for j in range (0,len(l2)):
        Sigma_color[j,:,:] = (F2[:,:,j]*F2[:,:,j].T)*sigma_cov[j]
    return Sigma_color
