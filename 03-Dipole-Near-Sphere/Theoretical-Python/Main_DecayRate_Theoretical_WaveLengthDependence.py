# %%
import matplotlib.pyplot as plt
import scipy.special as spc
import numpy as np


# The dispersion curve
def DrudeM(omega, omegaAu_D, gammaAu_D):
    varepsilon_D = omegaAu_D**2 / (omega**2 - 1j * (omega * gammaAu_D))
    return varepsilon_D


def LorentzM(omega, omegaAu_L, gammaAu_L, omega_0):
    varepsilon_L = omegaAu_L**2 / (omega**2 - omega_0**2 - 1j * omega * gammaAu_L)
    return varepsilon_L


def dispersion_Au(wavelength):  # , omegaAu_D, gammaAu_D, omegaAu_L, gammaAu_L, omega0_L):
    varepsilon_inf = 5.90157
    c_const = 3e8
    omegaAu_D = np.sqrt(1.297056e16**2 / varepsilon_inf)
    gammaAu_D = 4.108244e13
    omegaAu_L = np.sqrt(4.298305e15**2 * 1.26913 / varepsilon_inf)
    omega0_L = 4.298305e15
    # gammaAu_L = 2 * 4.104244e13 # Wrong Values Here!
    gammaAu_L = 0.374e15
    omega = 2 * np.pi * c_const / wavelength * 1e9
    varepsilon_gold = varepsilon_inf * \
        (1 - omegaAu_L**2 / (omega**2 - omega0_L**2 + 2 * 1j * omega * gammaAu_L) -
         omegaAu_D**2 / (omega**2 + 1j * (omega * gammaAu_D)))
    return varepsilon_gold

# Function to calculate the purcell factor


def sphereJ(l, x):
    return spc.spherical_jn(l, x)


def dphiJ(l, x):
    return spc.spherical_jn(l, x) * 0.5 + x * 0.5 * (spc.spherical_jn(l - 1, x) - spc.spherical_jn(l + 1, x))


def sphereH(l, x):
    l = l + 0.5
    return np.sqrt(np.pi / x / 2) * (spc.jv(l, x) + 1j * spc.yv(l, x))


def dphiH(l, x):
    return sphereH(l, x) * 0.5 + x * 0.5 * (sphereH(l - 1, x) - sphereH(l + 1, x))


def fun_nP_Cal(num_sum, kbr0, ksr0, index_s, kbr):
    rate_ortho = 0
    rate_para = 0
    for mm in range(num_sum):
        m = mm + 1
        a1 = (index_s**2 * sphereJ(m, ksr0) * dphiJ(m, kbr0) -
              sphereJ(m, kbr0) * dphiJ(m, ksr0))
        a2 = (index_s**2 * sphereJ(m, ksr0) * dphiH(m, kbr0) -
              sphereH(m, kbr0) * dphiJ(m, ksr0))
        b1 = (sphereJ(m, ksr0) * dphiJ(m, kbr0) -
              sphereJ(m, kbr0) * dphiJ(m, ksr0))
        b2 = (sphereJ(m, ksr0) * dphiH(m, kbr0) -
              sphereH(m, kbr0) * dphiJ(m, ksr0))
        bl = -a1 / a2
        al = -b1 / b2
        rate_ortho = rate_ortho + (2 * m + 1) * m * (m + 1) * \
            (bl) * (sphereH(m, kbr) / (kbr))**2
        rate_para = rate_para + (m + 0.5) * ((bl) * (dphiH(m, kbr) / (kbr))
                                             ** 2 + (al) * (sphereH(m, kbr))**2)
    rate_ortho = np.real(rate_ortho) * 3 / 2 + 1
    rate_para = 1 + 1.5 * np.real(rate_para)
    return rate_ortho, rate_para


# Define the struture and material

num = 100
lamda = np.linspace(300, 580, num) * 1e-9
index_b = 1
index_s = np.sqrt(dispersion_Au(lamda * 1e9))
k0 = 2 * np.pi / lamda
kb = k0 * index_b
ks = k0 * index_s

radius = 40e-9

# Calculate the theoretical expression
kbr0 = np.pi * 2 * radius / lamda
ksr0 = kbr0 * index_s
num_sum = 50
dis_theo = 10e-9 + radius
kbr = kb * dis_theo
nP_ortho_theo, nP_para_theo = fun_nP_Cal(
    num_sum, kbr0, ksr0, index_s, kbr)


fig1 = plt.figure()
plt.plot(lamda, nP_para_theo, 'r-', label="Tangential Dipole")
plt.plot(lamda, nP_ortho_theo, 'b', label="Radial Dipole")
plt.xlabel('distance-radius (nm)')
plt.ylabel('normalized power')
plt.title('theoretical Result')
plt.grid()
#plt.xlim([100e-9, 600e-9])
plt.legend()
# plt.savefig('../Figures/meep_verify.png')
plt.show()
