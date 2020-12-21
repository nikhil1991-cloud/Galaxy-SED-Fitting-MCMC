#Get indices for upper and lower data points near the predicted metallicity
metal_high = np.where(Metal>Given_Z)[0][0]
metal_low = np.where(Metal<Given_Z)[0][-1]
#Get indices for upper and lower data points near the predicted E(B-V)
EBV_I = Given_E/Delta_Ebv
EBV_high = int(np.floor(EBV_I))
EBV_low = int(np.floor(EBV_I+1))
#Interpolate for metallicity
SD4_1 = ((SD4[metal_low,:,:,:]*(Metal[metal_high]-Given_Z)) + (SD4[metal_high,:,:,:]*(Given_Z-Metal[metal_low])))/(Delta_metal)
SPR_1 = ((SPR[metal_low,:,:,:]*(Metal[metal_high]-Given_Z)) + (SPR[metal_high,:,:,:]*(Given_Z-Metal[metal_low])))/(Delta_metal)
SPB_1 = ((SPB[metal_low,:,:,:]*(Metal[metal_high]-Given_Z)) + (SPB[metal_high,:,:,:]*(Given_Z-Metal[metal_low])))/(Delta_metal)
SPI_1 = ((SPI[metal_low,:,:,:]*(Metal[metal_high]-Given_Z)) + (SPI[metal_high,:,:,:]*(Given_Z-Metal[metal_low])))/(Delta_metal)
SMg_1 = ((SMg[metal_low,2:,:,:]*(Metal[metal_high]-Given_Z)) + (SMg[metal_high,2:,:,:]*(Given_Z-Metal[metal_low])))/(Delta_metal)
SMgAB_1 = ((SMgAB[metal_low,2:,:]*(Metal[metal_high]-Given_Z)) + (SMgAB[metal_high,2:,:]*(Given_Z-Metal[metal_low])))/(Delta_metal)
#Interpolate for extinction
SD4_2 = ((SD4_1[:,EBV_low,:]*(EBV[EBV_high]-Given_E)) + (SD4_1[:,EBV_high,:]*(Given_E-EBV[EBV_low])))/(Delta_Ebv)
SPR_2 = (((SPR_1[:,EBV_low,:]*(EBV[EBV_high]-Given_E)) + (SPR_1[:,EBV_high,:]*(Given_E-EBV[EBV_low])))/(Delta_Ebv))
SPB_2 = (((SPB_1[:,EBV_low,:]*(EBV[EBV_high]-Given_E)) + (SPB_1[:,EBV_high,:]*(Given_E-EBV[EBV_low])))/(Delta_Ebv))
SPI_2 = (((SPI_1[:,EBV_low,:]*(EBV[EBV_high]-Given_E)) + (SPI_1[:,EBV_high,:]*(Given_E-EBV[EBV_low])))/(Delta_Ebv))
SMg_2 = (((SMg_1[:,EBV_low,:]*(EBV[EBV_high]-Given_E)) + (SMg_1[:,EBV_high,:]*(Given_E-EBV[EBV_low])))/(Delta_Ebv))
SMgAB_2 = ((SMgAB_1[:,EBV_low]*(EBV[EBV_high]-Given_E)) + (SMgAB_1[:,EBV_high]*(Given_E-EBV[EBV_low])))/(Delta_Ebv)
SMg_2[[0,1],:] = SMg_2[[1,0],:]
SMgAB_2[[0,1]] = SMgAB_2[[1,0]]

#For the interpolated SSP calculate the alpha matrices for regression
S_C0 = (SPR_2.T-SPB_2.T)*K_lick + SPB_2.T
Alpha_L= Lick_matrix(Lick_indices_data,S_C0,SPI_2.T,Delta_lmd)[Lick_list,:]
Alpha_D = D4000_matrix(Dn4000_data,SD4_2[1,:],SD4_2[0,:]).reshape(1,-1) #Dn4000
Alpha_C = Color_matrix(m1_list,m2_list,Mags_data,SMg_2.T,SMgAB_2,Ages)
All_alpha  = np.concatenate((Alpha_L,Alpha_D,Alpha_C),axis=0)
y = (-1)*All_alpha[:,-1]
A = All_alpha[:,0:(np.shape(All_alpha)[1]-1)]
alphas = np.zeros((9)) #initial guess
#generate sigma matrices for Lick, Dn4000 and Colors
Sigma_L = gen_lick_sigma(S_C0,Lick_indices_error,Lick_list,Ages)
Sigma_D = (gen_D4000_sigma(SD4_2[1,:],Dn4000_error,Ages)).reshape(1,10,10)
Sigma_C = gen_color_sigma(m1_list,m2_list,SMg_2,Mags_error,Mags_data,Ages)
ALL_sigma = np.concatenate((Sigma_L,Sigma_D,Sigma_C),axis=0)
#Join all sigmas to form sigma_double_dashed
Sigma_double_dashed = block_diag(*ALL_sigma)
#plt.imshow(Sigma_double_dashed,vmin=ALL_sigma.min(),vmax=ALL_sigma.max())
#Solve for alphas
alpha_bestfit = generate_alphas(A,y,alphas,Sigma_double_dashed,100)
#Creat SFH which is just an array of length 15 with 15 alphas
alpha_bestfit = np.clip(alpha_bestfit,0,alpha_bestfit.max())
alpha_bestfit = np.append(alpha_bestfit,1)
#Calculate Observables for CSP
D4000_Model = (np.sum(SD4_2[0,:]*alpha_bestfit))/(np.sum(SD4_2[1,:]*alpha_bestfit))
CR_T = np.sum(SPR_2*alpha_bestfit)
CB_T = np.sum(SPB_2*alpha_bestfit)
S0_T = (CR_T - CB_T)/(LR-LB)
C0_T = S0_T*(LI-LB) + CB_T
CI_T = np.sum(SPI_2*alpha_bestfit)
LICK_Model = (DI - CI_T/C0_T)
Photo = -2.5*np.log10(np.sum(SMg_2*alpha_bestfit)/(SMgAB_2))
MAGS_Model = Photo
