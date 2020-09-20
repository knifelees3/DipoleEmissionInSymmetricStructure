# BFP Image Program Function Needed


# The function list
"""
(1) Cal_RSP: To calculate the reflection coefficiencts for each layer
(2) Cal_Green: To calculate the Green function for a given layer and kx,ky
(3) Cal_Elec_Field: To calculate the electric field for a given p0 and Green
(4) Cal_Pattern: To calculate the pattern for a given electric field
(5) Cal_Dipole_From_Angle: To calculate the dipole p1,p2 from the angle (alpha,phi_1,phi_2)
(6) Grid_Data_TZH: To interploate the unstructured data from one coordinate to another
"""

# Import necessary package
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
# Basic Parameters
epsilon0 = 8.854187817620389850537e-12
mu0 = 1.2566370614359172953851e-6
const_c = np.sqrt(1 / epsilon0 / mu0)


# The Reflection and Transmission Function
# ______________________________________________________________________________________________________


def Cal_RSP(num_dl, Eplist, dl, nUp, nDn, klz):
    # The reflection coefficients
    RP21 = np.zeros(num_dl, dtype=complex)
    RP12 = np.zeros(num_dl, dtype=complex)
    RS21 = np.zeros(num_dl, dtype=complex)
    RS12 = np.zeros(num_dl, dtype=complex)

    for l in range(num_dl):
        RP21[l] = (Eplist[l] * klz[l + 1] - Eplist[l + 1] * klz[l]) / \
            (Eplist[l + 1] * klz[l] + Eplist[l] * klz[l + 1])
        RP12[l] = -RP21[l]
        RS21[l] = (klz[l + 1] - klz[l]) / (klz[l] + klz[l + 1])
        RS12[l] = -RS21[l]

    # For The Upper layer
    RPUp = 0 + 0 * 1j
    RSUp = 0 + 0 * 1j

    for m in range(nUp):
        l = num_dl - m - 1
        exp11 = np.exp(-1j * (klz[l + 1] - klz[l]) * dl[l])
        exp12 = np.exp(1j * (klz[l + 1] + klz[l]) * dl[l])
        exp21 = np.exp(-1j * (klz[l + 1] + klz[l]) * dl[l])
        exp22 = np.exp(1j * (klz[l + 1] - klz[l]) * dl[l])

        RPUp = (RPUp * exp11 + RP12[l] * exp12) / \
            (RPUp * RP12[l] * exp21 + exp22)
        RSUp = (RSUp * exp11 + RS12[l] * exp12) / \
            (RSUp * RS12[l] * exp21 + exp22)

    #  For the lower space
    RPDn = 0 + 0 * 1j
    RSDn = 0 + 0 * 1j
    for m in range(nDn):
        l = m + 1
        exp11 = np.exp(1j * (klz[l - 1] - klz[l]) * dl[l - 1])
        exp12 = np.exp(-1j * (klz[l - 1] + klz[l]) * dl[l - 1])
        exp21 = np.exp(+1j * (klz[l - 1] + klz[l]) * dl[l - 1])
        exp22 = np.exp(-1j * (klz[l - 1] - klz[l]) * dl[l - 1])

        RPDn = (RPDn * exp11 + RP21[l - 1] * exp12) / \
            (RPDn * RP21[l - 1] * exp21 + exp22)
        RSDn = (RSDn * exp11 + RS21[l - 1] * exp12) / \
            (RSDn * RS21[l - 1] * exp21 + exp22)

    return RSUp, RPUp, RSDn, RPDn, RS12, RP12, RS21, RP21


# The Green Function Will be Calculate By The Following Fucntion
# ______________________________________________________________________________________________________
def Cal_Green(nUp, nDn, num_layer, dl, POSD, dUpFar, dDnFar, RSUp, RPUp, RSDn, RPDn, RS12, RP12, RS21, RP21, kx, ky, klz):

    #  ************************************************************************************************************
    #  Initialize the matrix
    #  ************************************************************************************************************

    # Define the Green' tensor elements for each layer
    FSAUp = np.zeros((3, 3, nUp + 1), dtype=complex)
    FSBUp = np.zeros((3, 3, nUp + 1), dtype=complex)
    FPAUp = np.zeros((3, 3, nUp + 1), dtype=complex)
    FPBUp = np.zeros((3, 3, nUp + 1), dtype=complex)

    FSADn = np.zeros((3, 3, nDn + 1), dtype=complex)
    FSBDn = np.zeros((3, 3, nDn + 1), dtype=complex)
    FPADn = np.zeros((3, 3, nDn + 1), dtype=complex)
    FPBDn = np.zeros((3, 3, nDn + 1), dtype=complex)

    # Define the normalized tensor for s and p light
    SSUpA = np.zeros((3, 3, nUp + 1), dtype=complex)
    SSUpB = np.zeros((3, 3, nUp + 1), dtype=complex)
    PPUpA = np.zeros((3, 3, nUp + 1), dtype=complex)
    PPUpB = np.zeros((3, 3, nUp + 1), dtype=complex)

    SSDnA = np.zeros((3, 3, nDn + 1), dtype=complex)
    SSDnB = np.zeros((3, 3, nDn + 1), dtype=complex)
    PPDnA = np.zeros((3, 3, nDn + 1), dtype=complex)
    PPDnB = np.zeros((3, 3, nDn + 1), dtype=complex)

    # Define the coefficients for s and p light
    AUpS = np.zeros((3, 3, nUp + 1), dtype=complex)
    BUpS = np.zeros((3, 3, nUp + 1), dtype=complex)
    AUpP = np.zeros((3, 3, nUp + 1), dtype=complex)
    BUpP = np.zeros((3, 3, nUp + 1), dtype=complex)

    ADnS = np.zeros((3, 3, nDn + 1), dtype=complex)
    BDnS = np.zeros((3, 3, nDn + 1), dtype=complex)
    ADnP = np.zeros((3, 3, nDn + 1), dtype=complex)
    BDnP = np.zeros((3, 3, nDn + 1), dtype=complex)

    # Redifine the matrix for the upper layers and lower layers
    RS21Up = RS21[nDn:num_layer]
    RS12Up = RS12[nDn:num_layer]
    RP21Up = RP21[nDn:num_layer]
    RP12Up = RP12[nDn:num_layer]

    RS21Dn = RS21[0:nDn]
    RS12Dn = RS12[0:nDn]
    RP21Dn = RP21[0:nDn]
    RP12Dn = RP12[0:nDn]

    dlUp = np.vstack((POSD, dl[nDn:num_layer].reshape(nUp, 1)))
    dlDn = np.vstack((dl[0:nDn].reshape(nDn, 1), POSD))

    klzUp = klz[nDn:num_layer]
    klzDn = klz[0:nDn + 1]

    # The phase factor in the Green's tensor
    ZMatUp = dlUp
    ZMatDn = dlDn

    PhaseAUp = np.exp(1j * klzUp[:] * ZMatUp[:, 0] * 0)
    PhaseBUp = np.exp(-1j * klzUp[:] * ZMatUp[:, 0] * 0)

    PhaseADn = np.exp(1j * klzDn[:] * ZMatDn[:, 0] * 0)
    PhaseBDn = np.exp(-1j * klzDn[:] * ZMatDn[:, 0] * 0)
    #  ************************************************************************************************************
    #  Matrix for the s polarized light
    # *************************************************************************************************************

    # The matrix element is klz independent
    SSUpA[0, 0, :] = ky**2
    SSUpA[0, 1, :] = -kx * ky
    SSUpA[0, 2, :] = 0
    SSUpA[1, 0, :] = -kx * ky
    SSUpA[1, 1, :] = kx**2
    SSUpA[1, 2, :] = 0
    SSUpA[2, 0, :] = 0
    SSUpA[2, 1, :] = 0
    SSUpA[2, 2, :] = 0

    SSUpA = SSUpA / (kx**2 + ky**2)

    SSUpB[:, :, :] = SSUpA
    SSDnB[:, :, :] = SSUpA
    SSDnA[:, :, :] = SSUpA

    #  ************************************************************************************************************
    #  Matrix for the p polarized light
    # *************************************************************************************************************

    # The matrix elements is klz dependent
    # Upper layer
    for l in range(nUp + 1):
        PPUpA[0, 0, l] = kx**2 * klzUp[l]**2
        PPUpA[0, 1, l] = kx * ky * klzUp[l]**2
        PPUpA[0, 2, l] = -kx * klzUp[l] * (kx**2 + ky**2)
        PPUpA[1, 0, l] = kx * ky * klzUp[l]**2
        PPUpA[1, 1, l] = ky**2 * klzUp[l]**2
        PPUpA[1, 2, l] = -ky * klzUp[l] * (kx**2 + ky**2)
        PPUpA[2, 0, l] = -kx * klzUp[l] * (kx**2 + ky**2)
        PPUpA[2, 1, l] = -ky * klzUp[l] * (kx**2 + ky**2)
        PPUpA[2, 2, l] = (kx**2 + ky**2)**2
        PPUpA[:, :, l] = PPUpA[:, :, l] / \
            (kx**2 + ky**2 + klzUp[l]**2) / (kx**2 + ky**2)

    PPUpB = PPUpA
    PPUpB[0:2, :, :] = -PPUpB[0:2, :, :]

    for l in range(nDn + 1):
        PPDnA[0, 0, l] = -kx**2 * klzDn[l]**2
        PPDnA[0, 1, l] = -kx * ky * klzDn[l]**2
        PPDnA[0, 2, l] = +kx * klzDn[l] * (kx**2 + ky**2)
        PPDnA[1, 0, l] = -kx * ky * klzDn[l]**2
        PPDnA[1, 1, l] = -ky**2 * klzDn[l]**2
        PPDnA[1, 2, l] = +ky * klzDn[l] * (kx**2 + ky**2)
        PPDnA[2, 0, l] = -kx * klzDn[l] * (kx**2 + ky**2)
        PPDnA[2, 1, l] = -ky * klzDn[l] * (kx**2 + ky**2)
        PPDnA[2, 2, l] = (kx**2 + ky**2)**2
        PPDnA[:, :, l] = PPDnA[:, :, l] / \
            (kx**2 + ky**2 + klzDn[l]**2) / (kx**2 + ky**2)

    PPDnB = PPDnA
    PPDnB[0:2, :, :] = -PPDnB[0:2, :, :]

    #  ************************************************************************************************************
    #  Calculate the coefficients A_l,B_l in different layer
    # *************************************************************************************************************

    #  ************************************************************************************************************
    #  In the emitting layer
    # *************************************************************************************************************

    # For the s polarized light
    # #########The Upper Space############
    AUpS[0:2, 0:2, 0] = RSDn * (RSUp * np.exp(-1j * klzUp[0] * POSD) + np.exp(
        1j * klzUp[0] * POSD)) / (1 - RSUp * RSDn) + np.exp(-1j * klzUp[0] * POSD)

    BUpS[0:2, 0:2, 0] = RSUp * (RSDn * np.exp(1j * klzUp[0] * POSD) + np.exp(
        -1j * klzUp[0] * POSD)) / (1 - RSUp * RSDn)

    # #########The Lower Space#############

    ADnS[0:2, 0:2, nDn] = RSDn * (RSUp * np.exp(-1j * klzDn[nDn] * POSD) + np.exp(
        1j * klzDn[nDn] * POSD)) / (1 - RSUp * RSDn)

    BDnS[0:2, 0:2, nDn] = RSUp * (RSDn * np.exp(1j * klzDn[nDn] * POSD) + np.exp(
        -1j * klzDn[nDn] * POSD)) / (1 - RSUp * RSDn) + np.exp(1j * klzDn[nDn] * POSD)

    # For the P polarized light
    # #########The Upper Space#############
    AUpP[:, 0:2, 0] = RPDn * (RPUp * np.exp(-1j * klzUp[0] * POSD) - np.exp(
        1j * klzUp[0] * POSD)) / (1 - RPUp * RPDn) + np.exp(-1j * klzUp[0] * POSD)

    BUpP[:, 0:2, 0] = -RPUp * (RPDn * np.exp(1j * klzUp[0] * POSD) - np.exp(
        -1j * klzUp[0] * POSD)) / (1 - RPUp * RPDn)

    # For the third row, the sign and expressions should be carefully written,
    # since the sign is different
    AUpP[:, 2, 0] = RPDn * (RPUp * np.exp(-1j * klzUp[0] * POSD) + np.exp(
        1j * klzUp[0] * POSD)) / (1 - RPUp * RPDn) + np.exp(-1j * klzUp[0] * POSD)

    BUpP[:, 2, 0] = RPUp * (RPDn * np.exp(1j * klzUp[0] * POSD) + np.exp(
        -1j * klzUp[0] * POSD)) / (1 - RPUp * RPDn)

    # #########The Lower Space#############
    ADnP[:, 0:2, nDn] = -RPDn * (RPUp * np.exp(-1j * klzDn[nDn] * POSD) - np.exp(
        1j * klzDn[nDn] * POSD)) / (1 - RPUp * RPDn)

    BDnP[:, 0:2, nDn] = +RPUp * (RPDn * np.exp(1j * klzDn[nDn] * POSD) - np.exp(
        -1j * klzDn[nDn] * POSD)) / (1 - RPUp * RPDn) + np.exp(1j * klzDn[nDn] * POSD)

    # For the third row, the sign and expressions should be carefully written,
    # since the sign is different
    ADnP[:, 2, nDn] = RPDn * (RPUp * np.exp(-1j * klzDn[nDn] * POSD) + np.exp(
        1j * klzDn[nDn] * POSD)) / (1 - RPUp * RPDn)

    BDnP[:, 2, nDn] = RPUp * (RPDn * np.exp(1j * klzDn[nDn] * POSD) + np.exp(
        -1j * klzDn[nDn] * POSD)) / (1 - RPUp * RPDn) + np.exp(1j * klzDn[nDn] * POSD)

    #  ************************************************************************************************************
    #  In other layers
    # *************************************************************************************************************
    for l in range(nUp):

        expAA = np.exp(+1j * (klzUp[l] - klzUp[l + 1]) * dlUp[l + 1])
        expBA = np.exp(-1j * (klzUp[l] + klzUp[l + 1]) * dlUp[l + 1])
        expAB = np.exp(+1j * (klzUp[l] + klzUp[l + 1]) * dlUp[l + 1])
        expBB = np.exp(-1j * (klzUp[l] - klzUp[l + 1]) * dlUp[l + 1])

        AUpS[:, :, l + 1] = 1 / (1 - RS21Up[l]) * (AUpS[:, :, l]
                                                   * expAA + RS21Up[l] * BUpS[:, :, l] * expBA)
        BUpS[:, :, l + 1] = 1 / (1 - RS21Up[l]) * (RS21Up[l] * AUpS[:, :, l]
                                                   * expAB + BUpS[:, :, l] * expBB)

        AUpP[:, 0:2, l + 1] = 1 / (1 + RP21Up[l]) * (AUpP[:, 0:2, l]
                                                     * expAA + RP21Up[l] * BUpP[:, 0:2, l] * expBA)
        BUpP[:, 0:2, l + 1] = 1 / (1 + RP21Up[l]) * (RP21Up[l]
                                                     * AUpP[:, 0:2, l] * expAB + BUpP[:, 0:2, l] * expBB)

        AUpP[:, 2, l + 1] = klzUp[l + 1] / klzUp[l] / (1 + RP21Up[l]) * (
            AUpP[:, 2, l] * expAA + RP21Up[l] * BUpP[:, 2, l] * expBA)
        BUpP[:, 2, l + 1] = klzUp[l + 1] / klzUp[l] / (1 + RP21Up[l]) * (
            RP21Up[l] * AUpP[:, 2, l] * expAB + BUpP[:, 2, l] * expBB)

    for m in range(nDn):
        l = nDn - m
        expAA = np.exp(1j * (klzDn[l] - klzDn[l - 1]) * dlDn[l - 1])
        expBA = np.exp(-1j * (klzDn[l] + klzDn[l - 1]) * dlDn[l - 1])
        expAB = np.exp(1j * (klzDn[l] + klzDn[l - 1]) * dlDn[l - 1])
        expBB = np.exp(-1j * (klzDn[l] - klzDn[l - 1]) * dlDn[l - 1])

        ADnS[:, :, l - 1] = 1 / (1 - RS12Dn[l - 1]) * (ADnS[:, :, l]
                                                       * expAA + RS12Dn[l - 1] * BDnS[:, :, l] * expBA)
        BDnS[:, :, l - 1] = 1 / (1 - RS12Dn[l - 1]) * (RS12Dn[l - 1] * ADnS[:, :, l]
                                                       * expAB + BDnS[:, :, l] * expBB)

        ADnP[:, 0:2, l - 1] = 1 / (1 + RP12Dn[l - 1]) * (ADnP[:, 0:2, l]
                                                         * expAA + RP12Dn[l - 1] * BDnP[:, 0:2, l] * expBA)
        BDnP[:, 0:2, l - 1] = 1 / (1 + RP12Dn[l - 1]) * (RP12Dn[l - 1] * ADnP[:, 0:2, l]
                                                         * expAB + BDnP[:, 0:2, l] * expBB)

        ADnP[:, 2, l - 1] = klzDn[l - 1] / klzDn[l] / (1 + RP12Dn[l - 1]) * (
            ADnP[:, 2, l] * expAA + RP12Dn[l - 1] * BDnP[:, 2, l] * expBA)
        BDnP[:, 2, l - 1] = klzDn[l - 1] / klzDn[l] / (1 + RP12Dn[l - 1]) * (
            RP12Dn[l - 1] * ADnP[:, 2, l] * expAB + BDnP[:, 2, l] * expBB)

    #  ************************************************************************************************************
    #  Obtain the Green tensor
    # *************************************************************************************************************

    for l in range(nUp + 1):
        FSAUp[:, :, l] = SSUpA[:, :, l] * \
            AUpS[:, :, l] / klzUp[l] * PhaseAUp[l]
        FSBUp[:, :, l] = SSUpB[:, :, l] * \
            BUpS[:, :, l] / klzUp[l] * PhaseBUp[l]
        FPAUp[:, :, l] = PPUpA[:, :, l] * \
            AUpP[:, :, l] / klzUp[l] * PhaseAUp[l]
        FPBUp[:, :, l] = PPUpB[:, :, l] * \
            BUpP[:, :, l] / klzUp[l] * PhaseBUp[l]

    for m in range(nDn + 1):
        FSADn[:, :, m] = SSDnA[:, :, m] * \
            ADnS[:, :, m] / klzDn[m] * PhaseADn[m]
        FSBDn[:, :, m] = SSDnB[:, :, m] * \
            BDnS[:, :, m] / klzDn[m] * PhaseBDn[m]
        FPADn[:, :, m] = PPDnA[:, :, m] * \
            ADnP[:, :, m] / klzDn[m] * PhaseADn[m]
        FPBDn[:, :, m] = PPDnB[:, :, m] * \
            BDnP[:, :, m] / klzDn[m] * PhaseBDn[m]

    GreenSUp = (FSAUp + FSBUp) / 8 / np.pi**2 * 1j
    GreenPUp = (FPAUp + FPBUp) / 8 / np.pi**2 * 1j
    GreenSDn = (FSADn + FSBDn) / 8 / np.pi**2 * 1j
    GreenPDn = (FPADn + FPBDn) / 8 / np.pi**2 * 1j

    #  ************************************************************************************************************
    #  The far field green tensor
    # *************************************************************************************************************

    GreenSUp_Far = FSAUp[:, :, nUp] * \
        np.exp(1j * klzUp[nUp] * dUpFar) / 8 / np.pi**2 * 1j

    GreenPUp_Far = FPAUp[:, :, nUp] * \
        np.exp(1j * klzUp[nUp] * dUpFar) / 8 / np.pi**2 * 1j

    GreenSDn_Far = FSBDn[:, :, 0] * \
        np.exp(-1j * klzDn[0] * dDnFar) / 8 / np.pi**2 * 1j
    GreenPDn_Far = FPBDn[:, :, 0] * \
        np.exp(-1j * klzDn[0] * dDnFar) / 8 / np.pi**2 * 1j

    # return GreenSUp, GreenPUp, GreenSDn, GreenPDn, GreenSUp_Far,
    # GreenPUp_Far, GreenSDn_Far, GreenPDn_Far

    # Only return the far field part
    return GreenSUp_Far, GreenPUp_Far, GreenSDn_Far, GreenPDn_Far


# The electric field can be calculated by the following function
# For single green function
# ______________________________________________________________________________________________________
def Cal_Elec_Field(GreenSUp, GreenPUp, GreenSDn, GreenPDn, p0, omega, mu0):

    ESUp = 1j * omega * mu0 * np.matmul(GreenSUp, p0)
    EPUp = 1j * omega * mu0 * np.matmul(GreenPUp, p0)

    ESDn = 1j * omega * mu0 * np.matmul(GreenSDn, p0)
    EPDn = 1j * omega * mu0 * np.matmul(GreenPDn, p0)

    return ESUp, EPUp, ESDn, EPDn

# The pattern can be calculated by the following function
# ______________________________________________________________________________________________________


def Cal_Pattern(Eplist, epsilon0, mu0, num_layer, ESUp, EPUp, ESDn, EPDn, theta_Up, theta_Dn):

    CoeUp = 0.5 * np.real(np.sqrt(epsilon0 * Eplist[num_layer - 1] / mu0))
    CoeDn = 0.5 * np.real(np.sqrt(epsilon0 * Eplist[0] / mu0))

    PatternUpS = CoeUp * (np.abs(ESUp[0])**2 + np.abs(
        ESUp[1])**2 + np.abs(ESUp[2])**2) * np.abs(np.cos(theta_Up))**2

    PatternUpP = CoeUp * (np.abs(EPUp[0])**2 + np.abs(
        EPUp[1])**2 + np.abs(EPUp[2])**2) * np.abs(np.cos(theta_Up))**2

    PatternDnS = CoeDn * (np.abs(ESDn[0])**2 + np.abs(
        ESDn[1])**2 + np.abs(ESDn[2])**2) * np.abs(np.cos(theta_Dn))**2

    PatternDnP = CoeDn * (np.abs(EPDn[0])**2 + np.abs(
        EPDn[1])**2 + np.abs(EPDn[2])**2) * np.abs(np.cos(theta_Dn))**2

    PatternUp = PatternUpS + PatternUpP
    PatternDn = PatternDnS + PatternDnP

    return PatternUp, PatternDn,


# To calculate the two orthogonal dipole from angle
# ————————————————————————————————————————————————————————————————————————————————————————————————————————————————

def Cal_Dipole_From_Angle(angle):
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

#   To interpolate the original data into another type, I redefined the function so that
#   the input and output grid are all 2D which is more better for myself to use

# ________________________________________________________________________________________________


def Grid_Data_TZH(kx_grid, ky_grid, Values, kx_grid_in, ky_grid_in):
    kxy = np.transpose(np.vstack((np.ravel(kx_grid), np.ravel(ky_grid))))
    Values_1D = np.ravel(Values)
    Values_in = interpolate.griddata(kxy, Values_1D, (kx_grid_in, ky_grid_in), method='nearest')
    return Values_in

# _______________________________________________________________________________________________

# Plot Part
# to plot the single 2D map for


def Plot_Single_2D(kx_grid, ky_grid, Zmat, p0, name):

    fig_name = 'fig' + str(name)
    fig_name = plt.figure()
    plt.pcolormesh(kx_grid, ky_grid, abs(Zmat),
                   cmap='jet')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar()
    titlename = str(name) + ' ' + ' ' + 'p0=(' + \
        str(round(p0[0], 2)) + ',' + str(round(p0[1], 2)) + \
        ',' + str(round(p0[2], 2)) + ')'
    plt.title(titlename)
    filename = './Figures/picture_'  + '_' + str(name) + '_p0_' + \
        str(round(p0[0], 2)) + '_' + str(round(p0[1], 2)) + \
        '_' + str(round(p0[2], 2)) + '_.png'
    plt.savefig(filename)
    # plt.show()
    plt.close(fig_name)
    return 0

# to plot the single 1D map
# _______________________________________________________________________________________________


def Plot_Single_1D(kx, ymat, p0, name):
    fig_name = 'fig' + str(name)
    fig_name = plt.figure()
    plt.plot(kx, ymat, '-s')
    plt.xlabel('x')
    plt.ylabel('y')
    titlename = str(name) + ' ' + ' ' + 'p0=(' + \
        str(round(p0[0], 2)) + ',' + str(round(p0[1], 2)) + \
        ',' + str(round(p0[2], 2)) + ')'
    plt.title(titlename)
    filename = './Figures/picture_'  + '_' + str(name) + '_p0_' + \
        str(round(p0[0], 2)) + '_' + str(round(p0[1], 2)) + \
        '_' + str(round(p0[2], 2)) + '_.png'
    plt.savefig(filename)
    # plt.show()
    plt.close(fig_name)
    return 0

# To plot single 1D map and show the comparsion
# _______________________________________________________________________________________________


def Plot_Single_1D_Compare(kx, ymat_exp, ymat_fit, p0, name):
    fig_name = 'fig' + str(name)
    fig_name = plt.figure()
    plt.plot(kx, ymat_exp, 's')
    plt.plot(kx, ymat_fit, '-')
    plt.xlabel('x')
    plt.ylabel('y')
    titlename = str(name) + ' ' + ' ' + 'p0=(' + \
        str(round(p0[0], 2)) + ',' + str(round(p0[1], 2)) + \
        ',' + str(round(p0[2], 2)) + ')'
    plt.title(titlename)
    filename = './Figures/picture_'  + '_' + str(name) + '_p0_' + \
        str(round(p0[0], 2)) + '_' + str(round(p0[1], 2)) + \
        '_' + str(round(p0[2], 2)) + '_.png'
    plt.savefig(filename)
    # plt.show()
    plt.close(fig_name)
    return 0

# To plot four figures together
# _______________________________________________________________________________________________


def Plot_Whole_Compare(kx_grid, ky_grid, BFP_Exp, BFP_Fit, kx, ky, BFP_Exp_Rho, BFP_Exp_Phi, BFP_Fit_Rho, BFP_Fit_Phi, angle0):

    Plot_Single_2D(kx_grid, ky_grid, BFP_Exp, angle0, 'BFP_Exp')
    Plot_Single_2D(kx_grid, ky_grid, BFP_Fit, angle0, 'BFP_Fit')
    Plot_Single_1D_Compare(kx, BFP_Exp_Rho, BFP_Fit_Rho, angle0, 'BFP_Rho')
    Plot_Single_1D_Compare(ky, BFP_Exp_Phi, BFP_Fit_Phi, angle0, 'BFP_Phi')
    return 0

# To plot the structure we care


def Show_Structure(Eplist, dl, WL0, nUp, nDn, POSD):

    x = np.ones(10)

    # Obtain the thickness
    dis_ext = (dl[-1] - dl[0]) / (nUp + nDn + 1)
    dlext = np.ravel(np.vstack((dl[0] - dis_ext, dl[:], dl[-1] + dis_ext)))

    fig1 = plt.figure()
    for l in range(nUp + nDn):
        plt.plot(x * dl[l] * 1e9)
    plt.ylim((dlext[0] * 1e9, dlext[-1] * 1e9))

    for l in range(nUp + nDn + 1):
        strname = 'permitivity=' + str(Eplist[l, 0])
        plt.text(4, (dlext[l] + dlext[l + 1]) * 1e9 / 2, strname)
    plt.ylabel('Z (nm)')
    plt.xlabel('X Direction')
    plt.scatter(4, POSD * 1e9, s=5, marker="o")
    titlestr = 'Structure with dipole wavelength ' + str(WL0 * 1e9) + 'nm'
    plt.title(titlestr)
    filename = './Figures/Structure' + '.png'
    plt.savefig(filename)
    plt.show()

    return 0
