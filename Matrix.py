import numpy as np

'''Define Parameter matrix generating functions for Dn4000, Lick Indices, Colors'''

def Lick_matrix(Lick_array,Model_continuum,Model_index_flux,Delta_lmd):
    '''
    Input parameters:
          
          Lick array: (42,1) dimensional array containg all 42 lick indices from the data file
          
          Model_continuum: (j,42) dimensional array of 42 lick continuum of SSP models for j ages
          
          Model_index_flux: (j,42) dimensional array of 42 lick index flux of SSP models for j ages
    
          Delta_lmd: (42,1) diensional array of index band widths for lick indices
          
    Returns:
         
          (42,j) dimensional array of linear equations for lick indices
    '''
    return np.transpose(((Lick_array-Delta_lmd)*Model_continuum + Model_index_flux))
    
def D4000_matrix(D4000,Model_Blue,Model_Red):
    '''
    Input parameters:
    
        D4000: float64 Value of Dn4000 for the data
    
        Model_Blue: (j,1) dimensional array of blue band flux of models for j ages
    
        Model_Red: (j,1) dimensional array of red band flux of models for j ages
    
    Returns:
     
        (1,j) dimensional array of linear equation for Dn4000
    '''
    return np.transpose(((D4000*Model_Blue) - Model_Red))
    
def Color_matrix(l1,l2,magnitudes,models,models_AB,Ages):
    '''
    Input parameters:
        
        l1: List of indices of selected broad-band filter 1 (W1,W2,M2,u,g,r,i,z)
              
        l2: List of indices of selected broad-band filter 2 (W1,W2,M2,u,g,r,i,z)
              
        For example l1 = [0,1,2] and l2 = [3,3,3] will give 3 colors W1-u,W2-u and M2-u
        
        magnitudes: Array of shape (8,1) (W1,W2,M2,u,g,r,i,z) containing brodband magnitudes from data
              
        models: Array of shape (j,c) containing numerator of equation 7 from blanton et al for models at j ages
        
        models_AB: Array of shape (j,1) containing denominator of equation 7 from blanton et al for models at j ages
        
        Ages: Array of shape (1,j)
    
    Returns:
    
        (c,j) dimensional array of linear equations for colors
        
    '''
    color = magnitudes[l1] - magnitudes[l2]
    Color_matrix = models[:,l1]*(models_AB[l2]/models_AB[l1]) - models[:,l2]*(10**(-0.4*color))
    return Color_matrix.T
