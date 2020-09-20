import numpy as np
import matplotlib.pyplot as plt


WL0 = 580e-9
nUp = 2
nDn = 2

# Initialize the coordinate for each layer
# The dipole should be better in the coordinate z=0
dl = np.zeros((nUp + nDn, 1))
dis = 10e-9
dl[3] = 600e-9
dl[2] = 400e-9
dl[1] = 200e-9
dl[0] = 0

# The position of the dipole
POSD = dl[nDn - 1] + dis

# Dipole direction
# phi = np.pi / 2 + np.pi / 6
# dx = np.sin(16 / 180 * np.pi) * np.cos(phi)
# dy = np.sin(16 / 180 * np.pi) * np.sin(phi)
# dz = np.cos(16 / 180 * np.pi)
# p0 = np.array([dx, dy, dz])
p0 = np.array([0, 0, 1])

# Initialize the permitivity
Eplist = np.zeros((nUp + nDn + 1, 1), dtype=complex)

Eplist[4] = 1.78**2
Eplist[3] = 1.5**2
Eplist[2] = 1.78**2
Eplist[1] = 1.5**2
Eplist[0] = 1


def Show_Structure(Eplist, dl, WL0):

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
