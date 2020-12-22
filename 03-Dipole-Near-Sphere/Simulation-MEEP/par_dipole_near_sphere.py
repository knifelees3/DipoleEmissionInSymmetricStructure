"""
# Parallel version
# dipole near sphere orthogonal dipole (polarization perpendicular to radial axis)
# as a first step for learning python meep package
# date:2020 07 06
# author: Zhaohua Tian

"""
# %%
# from IPython import get_ipython

# %%
from matplotlib import pyplot as plt
import meep as mp
import numpy as np
# import PyMieScatt as ps
from meep.materials import Au
from tqdm import tqdm
import multiprocessing as mpr
import time


start = time.time()


# %%
# Define functions

# monitor for monitoring the out flux


def flux_monitor_set(sim, span, pos, frq):
    frq_cen = frq[0]
    dfrq = frq[1]
    nfrq = frq[2]

    x_span = mp.Vector3(0, span[1], span[2])
    y_span = mp.Vector3(span[0], 0, span[2])
    z_span = mp.Vector3(span[0], span[1], 0)

    pos_x1 = mp.Vector3(pos[0] - span[0] / 2, pos[1], pos[2])
    pos_x2 = mp.Vector3(pos[0] + span[0] / 2, pos[1], pos[2])
    pos_y1 = mp.Vector3(pos[0], pos[1] - span[1] / 2, pos[2])
    pos_y2 = mp.Vector3(pos[0], pos[1] + span[1] / 2, pos[2])
    pos_z1 = mp.Vector3(pos[0], pos[1], pos[2] - span[2] / 2)
    pos_z2 = mp.Vector3(pos[0], pos[1], pos[2] + span[2] / 2)

    box_x1 = sim.add_flux(frq_cen, dfrq, nfrq,
                          mp.FluxRegion(center=pos_x1, size=x_span))
    box_x2 = sim.add_flux(frq_cen, dfrq, nfrq,
                          mp.FluxRegion(center=pos_x2, size=x_span))
    box_y1 = sim.add_flux(frq_cen, dfrq, nfrq,
                          mp.FluxRegion(center=pos_y1, size=y_span))
    box_y2 = sim.add_flux(frq_cen, dfrq, nfrq,
                          mp.FluxRegion(center=pos_y2, size=y_span))
    box_z1 = sim.add_flux(frq_cen, dfrq, nfrq,
                          mp.FluxRegion(center=pos_z1, size=z_span))
    box_z2 = sim.add_flux(frq_cen, dfrq, nfrq,
                          mp.FluxRegion(center=pos_z2, size=z_span))
    return box_x1, box_x2, box_y1, box_y2, box_z1, box_z2

# %%
# monitor for sum up the power


def flux_monitor_cal(box_x1, box_x2, box_y1, box_y2, box_z1, box_z2):
    box_x1_flux = mp.get_fluxes(box_x1)
    box_x2_flux = mp.get_fluxes(box_x2)
    box_y1_flux = mp.get_fluxes(box_y1)
    box_y2_flux = mp.get_fluxes(box_y2)
    box_z1_flux = mp.get_fluxes(box_z1)
    box_z2_flux = mp.get_fluxes(box_z2)

    frqs = mp.get_flux_freqs(box_x1)
    p_flux = -(np.asarray(box_x1_flux) - np.asarray(box_x2_flux) + np.asarray(box_y1_flux) -
               np.asarray(box_y2_flux) + np.asarray(box_z1_flux) - np.asarray(box_z2_flux))
    return frqs, p_flux
# %%
# main calulations


def purcell_cal(sources, geometry,
                tstop, span, pos, frq):

    sim = mp.Simulation(cell_size=cell_size,
                        resolution=resolution,
                        boundary_layers=pml_layers,
                        sources=sources,
                        geometry=geometry)

    sim.run(until_after_sources=tstop)
    frqs, p_flux = flux_monitor_cal(
        box_x1, box_x2, box_y1, box_y2, box_z1, box_z2)
    return frqs, p_flux
# %%

# Define basci parameters
# _________________________________________________________________________


# Radius of sphere
r = 0.075  # um

# Define the pulse shape
wvl_min = 0.779  # um
wvl_max = 0.781  # um

frq_min = 1 / wvl_max
frq_max = 1 / wvl_min
frq_cen = 0.5 * (frq_min + frq_max)
dfrq = frq_max - frq_min
nfrq = 3
frq = [frq_cen, dfrq, nfrq]


# when to stop
tstop = 40

# Structure
resolution = 60  # pixels/μm
cell_size = mp.Vector3(3, 2, 2)
pml_layers = [mp.PML(thickness=0.5)]


# dipole polarize along x
pos = [0.0, 0.0, 0.0]
span = [r / 2, r / 2, r / 2]

sourcesx = [mp.Source(src=mp.GaussianSource(frq_cen, dfrq, nfrq),
                      center=mp.Vector3(pos[0], pos[1], pos[2]),
                      component=mp.Ex)]

# dipole polarize along z
sourcesz = [mp.Source(src=mp.GaussianSource(frq_cen, dfrq, nfrq),
                      center=mp.Vector3(pos[0], pos[1], pos[2]),
                      component=mp.Ez)]

# Define the sphere structure
dis_min = 0.1  # um
dis_max = 0.6  # um
num_dis = 4
dis_mat = np.linspace(dis_min, dis_max, num_dis)

# geometry0 = [mp.Sphere(material=mp.Medium(epsilon=1),
#                        center=mp.Vector3(x=r * 2),
#                        radius=r)]

# geometryAu = [mp.Sphere(material=Au,
#                         center=mp.Vector3(),
#                         radius=r)]

# %%
# Calculate the power 0
# print("Calculate the normalized power...")
# frqs, p0_flux = purcell_cal(cell_size, resolution, pml_layers,
#                             sourcesx, geometry0,
#                             tstop, span, pos, frq)

# np.savetxt("./data/p0_flux.txt", p0_flux)

# Calculate the power for different structure in a parallel way


def sweep_dis(l, sources, dis_mat,
              tstop, span, pos, frq):
    if sources == "x":
        sources_use = sourcesx
    if sources == "z":
        sources_use = sourcesz

    geometryAu = [mp.Sphere(material=Au,
                            center=mp.Vector3(x=r + dis_mat[l]),
                            radius=r)]
    frqs, p_flux = purcell_cal(sources_use, geometryAu,
                               tstop, span, pos, frq)
    return p_flux


if __name__ == '__main__':
    # for x dipole
    print("begin to calculate the x dipole")
    # pbar = tqdm(total=num_dis)
    pool = mpr.Pool(processes=2)
    results = [pool.apply_async(sweep_dis, args=(l, "x", dis_mat,
                                                 tstop, span, pos, frq
                                                 )) for l in range(num_dis)]
    pttmp = [p.get() for p in results]
    ptarray = np.asarray(pttmp)
    ptx = np.ravel(ptarray)
    # pool.terminate()
    pool.close()
    pool.join()

    # for z dipole
    print("begin to calculate the z dipole")
    # pbar = tqdm(total=num_dis)
    pool = mpr.Pool(processes=2)
    results = [pool.apply_async(sweep_dis, args=(l, "z", dis_mat,
                                                 tstop, span, pos, frq
                                                 )) for l in range(num_dis)]
    pttmp = [p.get() for p in results]
    ptarray = np.asarray(pttmp)
    ptz = np.ravel(ptarray)
    # pool.terminate()
    pool.close()
    pool.join()


end = time.time()
print("Total simulation time is:", str(end - start))

np.savetxt("./data/ptz_meep_ortho.txt", ptz)
np.savetxt("./data/ptx_meep_para.txt", ptx)
