import numpy as np
from numpy import linalg as LA
import sys

# my module
import molecule

# initiate class instances
m = molecule.Molecule()
x = molecule.Xray()
sa = molecule.Simulated_Annealing()

# command line arguments
xyzpool_file = str(sys.argv[1])  # name of xyzpool file (e.g. initconds.xyz)
nframes = int(sys.argv[2])  # number of frames to read from xyz trajectory
qtuple = (0.1, 8, 81)  # qmin, qmax, qlen
distance_indices = [0, 1, 2, 3, 4, 5]

def generate_pool(
    xyzpool_file,
    nframes,
    distance_indices,
    qtuple,
    elastic_bool=False,
    reference_xyz_file="reference.xyz",
):
    """
    creates: pool.npz
    """

    # read xyz trajectory
    print("reading xyz_traj file: %s" % xyzpool_file)
    natoms, comment, atomarray, xyzpool = m.read_xyz_traj(xyzpool_file, nframes)
    atomic_numbers = [m.periodic_table(symbol) for symbol in atomarray]

    nind = len(distance_indices)  # number of atoms distances will be calculated for
    dist_arr = np.zeros((nind, nind, nframes))
    for i in range(nframes):
        print("calculating distances for frame %i ..." % i)
        dist_arr[:, :, i] = m.distances_array(xyzpool[distance_indices, :, i])

    # definitions
    qmin, qmax, qlen = qtuple[0], qtuple[1], qtuple[2]
    qvector = np.linspace(qmin, qmax, qlen, endpoint=True)
    compton_array = x.compton_spline(atomic_numbers, qvector)

    # calculate reference IAM curve
    _, _, atomlist, reference_xyz = m.read_xyz(reference_xyz_file)
    reference_iam = x.iam_calc_compton(
        atomic_numbers, reference_xyz, qvector, compton_array, elastic_bool
    )

    # calculate pre-molecular IAM term
    aa, bb, cc = sa.read_iam_coeffs()
    compton, atomic_total, pre_molecular = sa.atomic_pre_molecular(
        atomic_numbers, qvector, aa, bb, cc
    )

    # calculate IAM and PCD for each frame
    pcd_arr = np.zeros((qlen, nframes))
    for i in range(nframes):
        print("calculating IAM and PCD for frame %i" % i)
        ##=#=#=# IAM CALCULATION #=#=#=##
        ## this takes 84% of the run time ... ##
        ## can it be optimised further? ##
        molecular = np.zeros(qlen)  # total molecular factor
        k = 0
        xyz_ = xyzpool[:, :, i]
        for ii in range(natoms):
            for jj in range(ii + 1, natoms):  # j > i
                qdij = qvector * LA.norm(xyz_[ii, :] - xyz_[jj, :])
                molecular += pre_molecular[k, :] * np.sin(qdij) / qdij
                k += 1
        iam_ = atomic_total + 2 * molecular
        if not elastic_bool:
            iam_ += compton
        ##=#=#=# END IAM CALCULATION #=#=#=##

        ##=#=#=# PCD & CHI2 CALCULATIONS #=#=#=##
        pcd_arr[:, i] = 100 * (iam_ / reference_iam - 1)


    #for i in range(nframes):
    #    print("calculating IAM and PCD for frame %i" % i)
    #    iam_ = x.iam_calc_compton(
    #        atomic_numbers, xyzpool[:, :, i], qvector, compton_array, elastic_bool
    #    )
    #    pcd_ = 100 * (iam_ / reference_iam - 1)
    #    # iam_arr[:, i] = iam_
    #    pcd_arr[:, i] = pcd_

    # Finally, save pool.npz
    np.savez(
        "pool.npz",
        natoms=natoms,
        atomarray=atomarray,
        atomic_numbers=atomic_numbers,
        nframes=nframes,
        xyzpool=xyzpool,
        distance_indices=distance_indices,
        dist_arr=dist_arr,
        reference_xyz=reference_xyz,
        reference_iam=reference_iam,
        qvector=qvector,
        elastic_bool=elastic_bool,
        pcd_arr=pcd_arr,
    )


# Run function
generate_pool(xyzpool_file, nframes, distance_indices, qtuple)

