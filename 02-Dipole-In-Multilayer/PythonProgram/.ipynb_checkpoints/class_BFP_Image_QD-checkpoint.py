
"""
# BFP Image Calculations of QD in Multi-layered structure

This is to calculate the BFP image for a given material and given
kx,ky distribution. Once kx,ky are given, then the pattern can be calclated.

# The description of this class can be found in the file "ReadMe.md"
"""

# %%
# Import commerical package
import numpy as np
import datetime
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Import my own package
import Fun_BFP_Image


# ______________________________________________________________________________
# %%
# The definition of the Class
class BFP_Image_QD:
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # **********************************Initialize Parameters*****************
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def __init__(self, Eplist, dl, nUp, nDn, p0, WL0, POSD):

        # Optical Constant
        self.epsilon0 = 8.854187817620389850537e-12
        self.mu0 = 1.2566370614359172953851e-6
        self.const_c = np.sqrt(1 / self.epsilon0 / self.mu0)

        # Main input parameters
        self.Eplist = Eplist
        self.dl = dl
        self.p0 = p0
        self.nUp = nUp
        self.nDn = nDn
        self.WL0 = WL0
        self.POSD = POSD

        # Test the input
        error_count = self.__check_input(Eplist, dl, nUp, nDn)
        if error_count == 0:
            # pass
            # f1 = open("process.txt", "a+")
            # f1.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') +
            #          ':The Basic Parameters Have Been Initialized!!! \n')
            # f1.close
            print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') +
                  ':The Basic Parameters Have Been Initialized!!!')
        else:
            # f1 = open("process.txt", "a+")
            # f1.write(datetime.datetime.now().strftime(
            #     '%Y-%m-%d %H:%M:%S') + ': WARNING: The parameters are not coorect, the results maybe not reliable!!! \n')
            # f1.close
            print(datetime.datetime.now().strftime(
                '%Y-%m-%d %H:%M:%S') + ': WARNING: The parameters are not coorect, the results maybe not reliable!!!')

        # Other parameters can be derived from the input parameters

        # number of layers
        self.num_layer = nUp + nDn + 1
        # number of coodinates
        self.num_dl = nUp + nDn
        # free space wave vector
        self.k0 = 2 * np.pi / WL0
        # Upper length to calculate the far field
        self.dUpFar = WL0 * 500
        # Lower length to calculate the far field
        self.dDnFar = -WL0 * 600
        # the wavevector in each layer
        self.kl = self.k0 * np.sqrt(Eplist)
        # the angular frequency
        self.omega = 2 * np.pi / self.WL0 * self.const_c

        # parameters that should be added laterly
        # The grid
        self.kx, self.ky, self.klz = 0, 0, 0
        # The size
        self.num_kx, self.num_ky = 0, 0
        # Green function initialize
        self.GreenSUp, self.GreenPUp, self.GreenSDn, self.GreenPDn = 0, 0, 0, 0
        # The angle list will be used
        self.theta_Up, self.theta_Dn = 0, 0

        # The clock number
        self.count = 0


# Check the correctness of the initialization

    def __check_input(self, Eplist, dl, nUp, nDn):
        size_Eplist = len(Eplist)
        size_dl = len(dl)

        error_count = 0
        if size_dl + 1 != size_Eplist:
            error_count = error_count + 1

        if nUp + nDn + 1 != size_Eplist:
            error_count = error_count + 1

        return error_count


# *******************************************************************
# ########################Functions Part#############################
# ___________________________________________________________________

# Part I Bacis Function Part
# *******************************************************************
# ___________________________________________________________________

# To calculate the R_{s/p} in different interface


    def Cal_RSP(self, klz):
        num_dl = self.num_dl
        Eplist = self.Eplist
        dl = self.dl
        nUp = self.nUp
        nDn = self.nDn

        RSUp, RPUp, RSDn, RPDn, RS12, RP12, RS21, RP21 = \
            Fun_BFP_Image.Cal_RSP(num_dl, Eplist, dl, nUp, nDn, klz)

        return RSUp, RPUp, RSDn, RPDn, RS12, RP12, RS21, RP21


# To calculate the Green function

    def Cal_Green(self, RSUp, RPUp, RSDn, RPDn, RS12, RP12, RS21, RP21, kx, ky, klz):
        nUp = self.nUp
        nDn = self.nDn
        num_layer = self.num_layer
        dl = self.dl
        POSD = self.POSD

        dUpFar = self.dUpFar
        dDnFar = self.dDnFar

        GreenSUp_Far, GreenPUp_Far, GreenSDn_Far, GreenPDn_Far = \
            Fun_BFP_Image.Cal_Green(nUp, nDn, num_layer, dl, POSD, dUpFar, dDnFar,
                                    RSUp, RPUp, RSDn, RPDn, RS12, RP12, RS21, RP21, kx, ky, klz)

        return GreenSUp_Far, GreenPUp_Far, GreenSDn_Far, GreenPDn_Far

# To calculate the electric field when the Green function is given
    def Cal_Elec_Field(self, GreenSUp, GreenPUp, GreenSDn, GreenPDn):
        p0 = self.p0
        omega = self.omega
        mu0 = self.mu0

        ESUp_Far, EPUp_Far, ESDn_Far, EPDn_Far = Fun_BFP_Image.Cal_Elec_Field(
            GreenSUp, GreenPUp, GreenSDn, GreenPDn,
            p0, omega, mu0)

        return ESUp_Far, EPUp_Far, ESDn_Far, EPDn_Far

# To calculate the emission pattern for a given elecric field
    def Cal_Pattern(self, ESUp, EPUp, ESDn, EPDn, theta_Up, theta_Dn):

        epsilon0 = self.epsilon0
        mu0 = self.mu0
        num_layer = self.num_layer
        Eplist = self.Eplist

        PatternUp, PatternDn = Fun_BFP_Image.Cal_Pattern(
            Eplist, epsilon0, mu0, num_layer, ESUp, EPUp, ESDn, EPDn, theta_Up, theta_Dn)

        nPatternUp = PatternUp / np.max(PatternUp)
        nPatternDn = PatternDn / np.max(PatternDn)
        return nPatternUp, nPatternDn

# To calculate the dipole orientation when the 2D dipole's angle is given
    def Cal_Dipole_From_Angle(self, angle):
        alpha = angle[0]
        phi_1 = angle[1]
        phi_2 = angle[2]
        # For Dipole 1
        d1x = np.sin(alpha) * np.cos(phi_1)
        d1y = np.sin(alpha) * np.sin(phi_1)
        d1z = np.cos(alpha)

        d2x = -np.cos(phi_2) * np.cos(alpha) * np.cos(phi_1) + \
            np.sin(phi_2) * np.sin(phi_1)

        d2y = -np.cos(phi_2) * np.cos(alpha) * np.sin(phi_1) - \
            np.sin(phi_2) * np.cos(phi_1)

        d2z = np.cos(phi_2) * np.sin(alpha)

        p1 = np.array([d1x, d1y, d1z])
        p2 = np.array([d2x, d2y, d2z])

        return p1, p2

# Part II Pattern Function Part
# ****************************************************************************
# __________________________________________________________________________________________
#   Single Pattern Cal with single kx,ky and a given p0
    def Cal_Pattern_Single_p0(self, kx, ky, p0):

        self.p0 = p0

        klz = np.sqrt(self.kl**2 - kx**2 - ky**2)

        theta_Up = np.arccos(klz[self.num_layer - 1] /
                             self.kl[self.num_layer - 1])
        theta_Dn = np.arccos(klz[0] / self.kl[0])

        RSUp, RPUp, RSDn, RPDn, RS12, RP12, RS21, RP21 = self.Cal_RSP(klz)

        GreenSUp, GreenPUp, GreenSDn, GreenPDn = self.Cal_Green(
            RSUp, RPUp, RSDn, RPDn, RS12, RP12, RS21, RP21, kx, ky, klz)

        ESUp, EPUp, ESDn, EPDn = self.Cal_Elec_Field(
            GreenSUp, GreenPUp, GreenSDn, GreenPDn)

        PatternUp, PatternDn = self.Cal_Pattern(
            ESUp, EPUp, ESDn, EPDn, theta_Up, theta_Dn)

        nPatternUp = PatternUp / np.max(PatternUp)
        nPatternDn = PatternDn / np.max(PatternDn)

        return nPatternUp, nPatternDn

#  Single Pattern Cal with single kx,ky and a given 2D angle

    def Cal_Pattern_single_QD(self, kx, ky, angle):
        p1, p2 = self.Cal_Dipole_From_Angle(angle)

        PatternUp1, PatternDn1 = self.Cal_Pattern_Single_p0(kx, ky, p1)
        PatternUp2, PatternDn2 = self.Cal_Pattern_Single_p0(kx, ky, p2)

        PatternUp = PatternUp1 + PatternUp2
        PatternDn = PatternDn1 + PatternDn2

        nPatternUp = PatternUp / np.max(PatternUp)
        nPatternDn = PatternDn / np.max(PatternDn)

        return PatternUp, PatternDn

    # _____________________________________________________________________________________________________
    # To calculate the green function for a given kx,ky list which can be used laterly

    def Cal_Green_List(self, kx, ky):
        self.kx, self.ky = kx, ky
        self.num_kx, self.num_ky = kx.shape

        klz = np.zeros((self.num_kx, self.num_ky,
                        self.num_layer), dtype=complex)
        for l in range(self.num_layer):
            klz[:, :, l] = np.sqrt(self.kl[l]**2 - self.kx**2 - self.ky**2)
        self.klz = klz

        self.theta_Up = np.arccos(self.klz[:, :, self.num_layer - 1] /
                                  self.kl[self.num_layer - 1])
        self.theta_Dn = np.arccos(self.klz[:, :, 0] / self.kl[0])

        GreenSUp, GreenPUp, GreenSDn, GreenPDn = np.zeros((3, 3, self.num_kx, self.num_ky), dtype=complex), np.zeros(
            (3, 3, self.num_kx, self.num_ky), dtype=complex), np.zeros((3, 3, self.num_kx, self.num_ky), dtype=complex), np.zeros((3, 3, self.num_kx, self.num_ky), dtype=complex)

        for l in range(self.num_kx):
            for m in range(self.num_ky):
                RSUp, RPUp, RSDn, RPDn, RS12, RP12, RS21, RP21 = self.Cal_RSP(
                    self.klz[l, m, :])
                GreenSUp[:, :, l, m], GreenPUp[:, :, l, m], GreenSDn[:, :, l, m], GreenPDn[:, :, l, m] = self.Cal_Green(
                    RSUp, RPUp, RSDn, RPDn, RS12, RP12, RS21, RP21, self.kx[l, m], self.ky[l, m], self.klz[l, m, :])
        self.GreenSUp, self.GreenPUp, self.GreenSDn, self.GreenPDn = GreenSUp, GreenPUp, GreenSDn, GreenPDn

        print(datetime.datetime.now().strftime(
            '%Y-%m-%d %H:%M:%S') + ': The Green Function Has Been Prepared')
        return 0

    # To calculate the pattern for a given dipole moment where the green function is prepared
    # For a given p1,p2
    def Cal_PatternUp_List_QD_p1p2(self, p1, p2):

        ESUp1 = np.zeros((self.num_kx, self.num_ky, 3), dtype=complex)
        EPUp1 = np.zeros((self.num_kx, self.num_ky, 3), dtype=complex)
        ESUp2 = np.zeros((self.num_kx, self.num_ky, 3), dtype=complex)
        EPUp2 = np.zeros((self.num_kx, self.num_ky, 3), dtype=complex)

        for l in range(3):
            ESUp1[:, :, l] = 1j * self.omega * self.mu0 * \
                (self.GreenSUp[l, 0, :, :] * p1[0] + self.GreenSUp[l,
                                                                   1, :, :] * p1[1] + self.GreenSUp[l, 2, :, :] * p1[2])
            EPUp1[:, :, l] = 1j * self.omega * self.mu0 * \
                (self.GreenPUp[l, 0, :, :] * p1[0] + self.GreenPUp[l,
                                                                   1, :, :] * p1[1] + self.GreenPUp[l, 2, :, :] * p1[2])
            ESUp2[:, :, l] = 1j * self.omega * self.mu0 * \
                (self.GreenSUp[l, 0, :, :] * p2[0] + self.GreenSUp[l,
                                                                   1, :, :] * p2[1] + self.GreenSUp[l, 2, :, :] * p2[2])
            EPUp2[:, :, l] = 1j * self.omega * self.mu0 * \
                (self.GreenPUp[l, 0, :, :] * p2[0] + self.GreenPUp[l,
                                                                   1, :, :] * p2[1] + self.GreenPUp[l, 2, :, :] * p2[2])

        # Cal the emission pattern
        PatternS = (np.abs(ESUp1[:, :, 0])**2 + np.abs(ESUp1[:, :, 1])**2 + np.abs(ESUp1[:, :, 2])**2 +
                    np.abs(ESUp2[:, :, 0])**2 + np.abs(ESUp2[:, :, 1])**2 + np.abs(ESUp2[:, :, 2])**2) * np.abs(np.cos(self.theta_Up[:, :]))**2
        PatternP = (np.abs(EPUp1[:, :, 0])**2 + np.abs(EPUp1[:, :, 1])**2 + np.abs(EPUp1[:, :, 2])**2 +
                    np.abs(EPUp2[:, :, 0])**2 + np.abs(EPUp2[:, :, 1])**2 + np.abs(EPUp2[:, :, 2])**2) * np.abs(np.cos(self.theta_Up[:, :]))**2
        Pattern = PatternS + PatternP
        nPattern = Pattern / np.max(Pattern)
        return nPattern



    # To calculate the pattern for a given dipole moment where the green function is prepared
    # For a given p1
    # The upper and lower pattern will both be simulated
    def Cal_Pattern_List_QD_p1(self, p1):

        ESUp = np.zeros((self.num_kx, self.num_ky, 3), dtype=complex)
        EPUp = np.zeros((self.num_kx, self.num_ky, 3), dtype=complex)
        ESDn = np.zeros((self.num_kx, self.num_ky, 3), dtype=complex)
        EPDn = np.zeros((self.num_kx, self.num_ky, 3), dtype=complex)

        for l in range(3):
            ESUp[:, :, l] = 1j * self.omega * self.mu0 * \
                (self.GreenSUp[l, 0, :, :] * p1[0] + self.GreenSUp[l,
                                                                   1, :, :] * p1[1] + self.GreenSUp[l, 2, :, :] * p1[2])
            EPUp[:, :, l] = 1j * self.omega * self.mu0 * \
                (self.GreenPUp[l, 0, :, :] * p1[0] + self.GreenPUp[l,
                                                                   1, :, :] * p1[1] + self.GreenPUp[l, 2, :, :] * p1[2])
            ESDn[:, :, l] = 1j * self.omega * self.mu0 * \
                (self.GreenSDn[l, 0, :, :] * p1[0] + self.GreenSDn[l,
                                                                   1, :, :] * p1[1] + self.GreenSDn[l, 2, :, :] * p1[2])
            EPDn[:, :, l] = 1j * self.omega * self.mu0 * \
                (self.GreenPDn[l, 0, :, :] * p1[0] + self.GreenPDn[l,
                                                                   1, :, :] * p1[1] + self.GreenPDn[l, 2, :, :] * p1[2])


        # Cal the emission pattern
        PatternUpS = (np.abs(ESUp[:, :, 0])**2 + np.abs(ESUp[:, :, 1])**2 + np.abs(ESUp[:, :, 2])**2) * np.abs(np.cos(self.theta_Up[:, :]))**2
        PatternUpP = (np.abs(EPUp[:, :, 0])**2 + np.abs(EPUp[:, :, 1])**2 + np.abs(EPUp[:, :, 2])**2) * np.abs(np.cos(self.theta_Up[:, :]))**2
        PatternDnS = (np.abs(ESDn[:, :, 0])**2 + np.abs(ESDn[:, :, 1])**2 + np.abs(ESDn[:, :, 2])**2) * np.abs(np.cos(self.theta_Dn[:, :]))**2
        PatternDnP = (np.abs(EPDn[:, :, 0])**2 + np.abs(EPDn[:, :, 1])**2 + np.abs(EPDn[:, :, 2])**2) * np.abs(np.cos(self.theta_Dn[:, :]))**2
        PatternUp = PatternUpS + PatternUpP
        PatternDn = PatternDnS + PatternDnP
        nPatternUp = PatternUp 
        nPatternDn = PatternDn
        return nPatternUp,nPatternDn

    # For a given angle
    def Cal_PatternUp_List_QD_Angle(self, alpha, phi_1, phi_2):
        angle = [alpha, phi_1, phi_2]
        p1, p2 = self.Cal_Dipole_From_Angle(angle)
        Pattern = self.Cal_PatternUp_List_QD_p1p2(p1, p2)
        nPattern = Pattern / np.max(Pattern)
        return nPattern

        # Part III Fitting Function Part
        # ***********************************************************************
# ____________________________________________________________________________________________________
    # I use alpha,phi_1,phi_2 here because the parameters should be single rather than array

    # The function used to fit , the function will use "Cal_Pattern_Single_p0" to calculate which
    # is less efficienct
    def Cal_Pattern_Single_QD_Fit(self, kxy, alpha, phi_1, phi_2):
        angle = [alpha, phi_1, phi_2]
        numk, num_temp = kxy.shape
        p1, p2 = self.Cal_Dipole_From_Angle(angle)

        PatternUp, PatternDn = np.zeros(numk), np.zeros(numk)

        self.count = self.count + 1
        print(datetime.datetime.now().strftime(
            '%Y-%m-%d %H:%M:%S') + ' Iteration Step ', self.count)

        for l in range(numk):

            kx, ky = kxy[l, 0], kxy[l, 1]

            PatternUp1, PatternDn1 = self.Cal_Pattern_Single_p0(kx, ky, p1)
            PatternUp2, PatternDn2 = self.Cal_Pattern_Single_p0(kx, ky, p2)

            PatternUp[l] = PatternUp1 + PatternUp2
            PatternDn[l] = PatternDn1 + PatternDn2

        # Normalize the data
        nPUp = np.abs(PatternUp) / np.max(np.abs(PatternUp))
        # nPDn = np.abs(PatternDn) / np.max(np.abs(PatternDn))

        # I will use the down direction as the BFP image direction
        nPUp1d = np.ravel(nPUp)
        # nPDn1d = np.ravel(nPDn)
        return nPUp1d

    # Fit Function (Less Efficient)
    def BFP_Curvefit_Single_Angle(self, kxy, Exp_data, angle0, para_bounds):
        self.count = 0
        print(datetime.datetime.now().strftime(
            '%Y-%m-%d %H:%M:%S') + ': Begin To Fit')
        # perform the fit, making sure to flatten the noisy data for the fit routine
        fit_params, cov_mat = curve_fit(self.Cal_Pattern_Single_QD_Fit,
                                        kxy, Exp_data, p0=angle0, bounds=para_bounds)

        # calculate fit parameter errors from covariance matrix
        fit_errors = np.sqrt(np.diag(cov_mat))

        # manually calculate R-squared goodness of fit
        fit_residual = Exp_data - \
            self.Cal_Pattern_Single_QD_Fit(kxy, *fit_params)
        fit_Rsquared = 1 - np.var(fit_residual) / np.var(Exp_data)

        print('Fit R-squared:', fit_Rsquared, '\n')
        print('Fit $\alpha$ :', fit_params[0], '\u00b1', fit_errors[0])
        print('Fit $\phi_1$: ', fit_params[1], '\u00b1', fit_errors[1])
        print('Fit $\phi_2$: ', fit_params[2], '\u00b1', fit_errors[2])

        return fit_params, fit_errors

    def Cal_Pattern_List_QD_Fit_Angle(self, kxy, alpha, phi_1, phi_2):
        # Fit Show
        self.count = self.count + 1
        print(datetime.datetime.now().strftime(
            '%Y-%m-%d %H:%M:%S') + ': Iteration Step ', self.count)

        # We will ignore kxy actually
        Pattern = self.Cal_PatternUp_List_QD_Angle(alpha, phi_2, phi_1)
        nPattern = Pattern / np.max(Pattern)
        nPattern1D = np.ravel(nPattern)

        return nPattern1D

    def BFP_Curvefit_List_Angle(self, kxy, Exp_data, angle0, para_bounds):

        print(datetime.datetime.now().strftime(
            '%Y-%m-%d %H:%M:%S') + ': Begin To Fit')
        # perform the fit, making sure to flatten the noisy data for the fit routine
        self.count = 0
        fit_params, cov_mat = curve_fit(self.Cal_Pattern_List_QD_Fit_Angle,
                                        kxy, Exp_data, p0=angle0, bounds=para_bounds)

        # calculate fit parameter errors from covariance matrix
        fit_errors = np.sqrt(np.diag(cov_mat))
        # manually calculate R-squared goodness of fit
        fit_residual = Exp_data - \
            self.Cal_Pattern_List_QD_Fit_Angle(kxy, *fit_params)
        fit_Rsquared = 1 - np.var(fit_residual) / np.var(Exp_data)

        print('Fit R-squared:', fit_Rsquared, '\n')
        print('Fit alpha :', fit_params[0], '\u00b1', fit_errors[0])
        print('Fit phi_1: ', fit_params[1], '\u00b1', fit_errors[1])
        print('Fit phi_2: ', fit_params[2], '\u00b1', fit_errors[2])

        return fit_params, fit_errors

    # ______________________________________________________________________________________________

    # Part IV Data Process Part
    # **********************************************************************************************
    def Cal_RhoPhi_Dis(self, Pattern_RhoPhi):
        Pattern_Rho = np.sum(Pattern_RhoPhi, axis=0)
        Pattern_Phi = np.sum(Pattern_RhoPhi, axis=1)
        return Pattern_Rho, Pattern_Phi

    def Trans_XY_to_RhoPhi(self, kx_grid, ky_grid, Pattern_XY, kx_grid_in, ky_grid_in):
        Pattern_RhoPhi = Fun_BFP_Image.Grid_Data_TZH(
            kx_grid, ky_grid, Pattern_XY, kx_grid_in, ky_grid_in)
        Pattern_Rho, Pattern_Phi = self.Cal_RhoPhi_Dis(Pattern_RhoPhi)

        return Pattern_Rho, Pattern_Phi

    # Part V Structure visualization part
    def Show_Structure(self):
        Fun_BFP_Image.Show_Structure(self.Eplist, self.dl, self.WL0, self.nUp, self.nDn, self.POSD)
        
    # Add in 2021 07 27
    # The definition of Dipole3D and 3D pattern
    def Dipole3D(self,para):
        alpha=para[0]
        beta=para[1]
        phix=para[2]
        phiy=para[3]
        phiz=para[4]
        dx=np.transpose([1,0,0])
        dy=np.transpose([0,alpha,0])
        dz=np.transpose([0,0,beta])

        norm=np.sqrt(1**2+alpha**2+beta**2)

        rotation=np.zeros((3,3))
        rotation[:,0]=np.transpose(np.array([np.cos(phiz)*np.cos(phiy),np.sin(phiz)*np.cos(phiy),-np.sin(phiy)]))
        rotation[:,1]=np.transpose(np.array([np.cos(phiz)*np.sin(phiy)*np.sin(phix)-np.sin(phiz)*np.cos(phix),np.sin(phiz)*np.sin(phiy)*np.sin(phix)+np.cos(phiz)*np.cos(phix),np.cos(phiy)*np.sin(phix)]))
        rotation[:,2]=np.transpose(np.array([np.cos(phiz)*np.sin(phiy)*np.cos(phix)+np.sin(phiz)*np.sin(phix),np.sin(phiz)*np.sin(phiy)*np.cos(phix)-np.cos(phiz)*np.sin(phix),np.cos(phiy)*np.cos(phix)]))

        d1=np.dot(rotation,dx)/norm
        d2=np.dot(rotation,dy)/norm
        d3=np.dot(rotation,dz)/norm

        return d1,d2,d3

    def Pattern3D(self,para):
        d1,d2,d3=self.Dipole3D(para)
        PatternUpd1,PatternDnd1 = self.Cal_Pattern_List_QD_p1(d1)
        PatternUpd2,PatternDnd2 = self.Cal_Pattern_List_QD_p1(d2)
        PatternUpd3,PatternDnd3 = self.Cal_Pattern_List_QD_p1(d3)
        PatternUpd1[np.isnan(PatternUpd1)]=0
        PatternUpd2[np.isnan(PatternUpd2)]=0
        PatternUpd3[np.isnan(PatternUpd3)]=0
        PatternDnd1[np.isnan(PatternDnd1)]=0
        PatternDnd2[np.isnan(PatternDnd2)]=0
        PatternDnd3[np.isnan(PatternDnd3)]=0
        PatternUp=PatternUpd1+PatternUpd2+PatternUpd3
        PatternDn=PatternDnd1+PatternDnd2+PatternDnd3
        return PatternUp,PatternDn
