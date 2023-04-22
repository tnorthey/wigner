import numpy as np
import sys
from timeit import default_timer
# my module
import molecule

start = default_timer()

# command line arguments
xyztraj_file = str(sys.argv[1]) # name of xyztraj file
nframes = int(sys.argv[2])      # number of frames to read from xyz trajectory
distances = bool(sys.argv[3])   # bool, calculate all the distances or not
xray = bool(sys.argv[4])        # bool, calculate all percent diff or not

# initiate class instances
m = molecule.Molecule()
x = molecule.Xray()

# definitions

# read xyz trajectory
natoms, comment, atomarray, xyztraj = m.read_xyz_traj(xyztraj_file, nframes)
atomic_numbers = [m.periodic_table(symbol) for symbol in atomarray]
print(atomarray)

if distances:
    # distances
    indices_to_save = [0, 1, 2, 3, 4, 5]  # these are the carbons in CHD 
    print('Remember: modify indices if changing molecule!')

    dist_arr = np.zeros((natoms, natoms, nframes))
    for i in range(nframes):
        dist_arr[:, :, i] = m.distances_array(xyztraj[:, :, i])
    
    # save to dat files
    for i in range(len(indices_to_save)):
        for j in range(i + 1, len(indices_to_save)):
            np.savetxt('distances/r%i%i.dat' % (i + 1, j + 1), dist_arr[i, j, :] )

if xray:
    
    # definitions
    qmin = 0.1
    qmax = 8.0
    qlen = 81
    qvector = np.linspace(qmin, qmax, qlen, endpoint=True)
    elastic = False
    compton_array = x.compton_spline(atomic_numbers, qvector)
    
    # calculate reference IAM curve
    reference_xyz_file = 'xyz/equilibrium.xyz'
    _, _, atomlist, reference_xyz = m.read_xyz(reference_xyz_file)
    reference_iam = x.iam_calc_compton(atomic_numbers, reference_xyz, qvector, compton_array, elastic)

    # read target PCD file
    # for now just create it from the xyz.
    # In general, I should read in pure data (the xyz wouldn't be known)
    target_xyz_file = 'xyz/target.xyz'
    _, _, atomlist, target_xyz = m.read_xyz(target_xyz_file)
    target_iam = x.iam_calc_compton(atomic_numbers, target_xyz, qvector, compton_array, elastic)
    target_pcd = 100 * (target_iam / reference_iam - 1)
    # save to file, why not
    np.savetxt('pcd/target.dat', np.column_stack((qvector, target_pcd)) )
    
    # calculate IAM and PCD for each frame
    iam_arr = np.zeros((qlen, nframes))
    pcd_arr = np.zeros((qlen, nframes))
    chi2_arr = np.zeros(nframes)
    for i in range(nframes):
        iam_ = x.iam_calc_compton(atomic_numbers, xyztraj[:, :, i], qvector, compton_array, elastic)
        iam_arr[:, i] = iam_
        pcd_ = 100 * (iam_ / reference_iam - 1)
        pcd_arr[:, i] = pcd_
        # chi2
        chi2_arr[i] = np.sum((pcd_ - target_pcd) ** 2) / qlen
        np.savetxt('pcd/%s_sample.dat' % str(i).zfill(6), np.column_stack((qvector, pcd_)) )

# save chi2_array
np.savetxt('pcd/chi2_sample.dat', chi2_arr)

print( "Total time: %3.2f s" % float(default_timer() - start) )

