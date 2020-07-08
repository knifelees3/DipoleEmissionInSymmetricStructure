# %%
import matplotlib.pyplot as plt
import scipy.special as spc
import numpy as np


def sphereJ(l, x):
    return spc.spherical_jn(l, x)


def dphiJ(l, x):
    return spc.spherical_jn(l, x) * 0.5 + x * 0.5 * (spc.spherical_jn(l - 1, x) - spc.spherical_jn(l + 1, x))


def sphereH(l, x):
    l = l + 0.5
    return np.sqrt(np.pi / x / 2) * (spc.jv(l, x) + 1j * spc.yv(l, x))


def dphiH(l, x):
    return sphereH(l, x) * 0.5 + x * 0.5 * (sphereH(l - 1, x) - sphereH(l + 1, x))


def fun_nP_Cal(num_sum, num_dis, dis, alpha, alpham, index_s, k):
    nP_para = np.zeros(num_dis)
    nP_ortho = np.zeros(num_dis)
    for l in range(num_dis):
        rate_ortho = 0
        rate_para = 0
        r = dis[l]
        for mm in range(num_sum):
            m = mm + 1
            a1 = (index_s**2 * sphereJ(m, alpham) * dphiJ(m, alpha) -
                  sphereJ(m, alpha) * dphiJ(m, alpham))
            a2 = (index_s**2 * sphereJ(m, alpham) * dphiH(m, alpha) -
                  sphereH(m, alpha) * dphiJ(m, alpham))
            b1 = (sphereJ(m, alpham) * dphiJ(m, alpha) -
                  sphereJ(m, alpha) * dphiJ(m, alpham))
            b2 = (sphereJ(m, alpham) * dphiH(m, alpha) -
                  sphereH(m, alpha) * dphiJ(m, alpham))
            x = k * r
            al = a1 / a2
            bl = b1 / b2
            rate_ortho = rate_ortho + (2 * m + 1) * m * (m + 1) * \
                (-al) * (sphereH(m, x) / (x))**2
            rate_para = rate_para + (m + 0.5) * ((-al) * (dphiH(m, x) / (x))
                                                 ** 2 + (-bl) * (sphereH(m, x))**2)
        rate_ortho = np.real(rate_ortho) * 3 / 2 + 1
        rate_para = 1 + 1.5 * np.real(rate_para)
        nP_para[l] = rate_para
        nP_ortho[l] = rate_ortho
    return nP_para, nP_ortho


# Define the struture and material
index_s = np.sqrt(-22.473 - 1.3974 * 1j)
lamda = 780e-9
k0 = 2 * np.pi / lamda
k = k0
radius = 75e-9

# Calculate the theoretical expression
alpha = np.pi * 2 * radius / lamda
alpham = alpha * index_s
num_dis = 100
num_sum = 40
dis_theo = np.linspace(radius + 50e-9, radius + 600e-9, num_dis)
nP_ortho_theo, nP_para_theo = fun_nP_Cal(
    num_sum, num_dis, dis_theo, alpha, alpham, index_s, k)
# %%
# Load the data from meep
p0_data = np.loadtxt("../Data/p0_meep.txt")
ptx_data = np.loadtxt("../Data/ptx_meep_para_new.txt")
ptz_data = np.loadtxt("../Data/ptz_meep_ortho_new.txt")
dis_meep = np.linspace(50e-9, 600e-9, 22) + radius
ptpara_mat = np.reshape(ptx_data, (22, 3))
ptortho_mat = np.reshape(ptz_data, (22, 3))

nP_para_meep = np.zeros((22, 3))
nP_ortho_meep = np.zeros((22, 3))
for l in range(3):
    nP_para_meep[:, l] = ptpara_mat[:, l] / p0_data[l]
    nP_ortho_meep[:, l] = ptortho_mat[:, l] / p0_data[l]
# %%
import matplotlib.pyplot as plt
fig1 = plt.figure()
plt.plot(dis_theo - radius, nP_para_theo, 'r-', label="Parallel Dipole Theory")
plt.plot(dis_theo - radius, nP_ortho_theo, 'b', label="Orthogonal Dipole Theory")
plt.plot(dis_meep - radius, nP_para_meep[:, 1], 'r>', label="Parallel Dipole Meep")
plt.plot(dis_meep - radius, nP_ortho_meep[:, 1], 'bo', label="Orthogonal Dipole Meep")
plt.xlabel('distance-radius (nm)')
plt.ylabel('normalized power')
plt.title('MEEP Simulation Result')
plt.grid()
plt.xlim([100e-9, 600e-9])
plt.legend()
plt.savefig('../Figures/meep_verify.png')
plt.show()


# %%
