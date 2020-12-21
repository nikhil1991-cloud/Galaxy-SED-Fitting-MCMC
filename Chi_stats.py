import numpy as np

"""These functions are used to calculate chi-square between data and models observables"""


def chi_square(Data,Model,Data_sigma):
    '''
    This function is used to calculate the chi-square for Dn4000 and Lick indices.
    '''
    return ((Data - Model)/(Data_sigma))**2
    
def chi_square_color(l1,l2,Data_mag,Model_mag,Data_mag_sigma):
    '''
    This function is used to calculate chi-square for c colors.
    
    Input Parameters:
          
          l1: List of indices of selected broad-band filter 1 (W1,W2,M2,u,g,r,i,z).
                
          l2: List of indices of selected broad-band filter 2 (W1,W2,M2,u,g,r,i,z).
                
          For example l1 = [0,1,2] and l2 = [3,3,3] will give 3 colors W1-u,W2-u and M2-u.
          
          Data_mag: Array of shape (8,1) (W1,W2,M2,u,g,r,i,z) containing brodband magnitudes from data.
                
          Model_mag: Array of shape (8,1) containing model 8 magnitudes.
          
          Data_mag_sigma: Array of shape (8,1) (W1,W2,M2,u,g,r,i,z) containing 1 sigma values of brodband magnitudes from data.
    '''

    chi_sq_color = np.zeros(len(l1))
    Color_error = (np.sqrt(Data_mag_sigma[l1]**2 + Data_mag_sigma[l2]**2))
    Model_color = Model_mag[l1] - Model_mag[l2]
    Data_color = Data_mag[l1] - Data_mag[l2]
    chi_sq_color = ((Data_color - Model_color)/Color_error)**2
    return chi_sq_color
