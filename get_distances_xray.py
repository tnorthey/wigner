import numpy as np
import sys
from timeit import default_timer

# my module
import molecule

start = default_timer()

# command line arguments
xyztraj_file = str(sys.argv[1])  # name of xyztraj file
nframes = int(sys.argv[2])  # number of frames to read from xyz trajectory
distances = bool(sys.argv[3])  # bool, calculate all the distances or not
xray = bool(sys.argv[4])  # bool, calculate all percent diff or not

# initiate class instances
m = molecule.Molecule()
x = molecule.Xray()

# definitions
save_individual_pcd = False

# read xyz trajectory
## add: check if file exists
xyztraj_npz = True
if xyztraj_npz:
    data = np.load("xyztraj.npz")
    natoms = data["natoms"]
    atomarray = data["atomarray"]
    atomic_numbers = data["atomic_numbers"]
    nframes = data["nframes"]
    xyztraj = data["xyztraj"]
else:
    print('reading xyz_traj file: %s' % xyztraj_file)
    natoms, comment, atomarray, xyztraj = m.read_xyz_traj(xyztraj_file, nframes)
    atomic_numbers = [m.periodic_table(symbol) for symbol in atomarray]
    np.savez(
        "xyztraj.npz",
        natoms=natoms,
        atomarray=atomarray,
        atomic_numbers=atomic_numbers,
        nframes=nframes,
        xyztraj=xyztraj,
    )

if distances:
    # distances
    indices_to_save = [0, 1, 2, 3, 4, 5]  # these are the carbons in CHD
    print("Remember: modify indices if changing molecule!")

    nind = len(indices_to_save)
    dist_arr = np.zeros((nind, nind, nframes))
    for i in range(nframes):
        ## is this slow? It might be fine to leave it if it's fast enough.
        ## add: only calculate indices specified, not all n(n-1)/2 distances...
        print("calculating distances for frame %i ..." % i)
        #print(xyztraj[indices_to_save, :, i])
        dist_arr[:, :, i] = m.distances_array(
            xyztraj[indices_to_save, :, i]
        )

    # save to dat files
    print("saving distances to dat files")
    for i in range(nind):
        for j in range(i + 1, nind):
            print('saving distance %i %i to file...' % (i, j))
            np.savetxt("distances/r%i%i.dat" % (i + 1, j + 1), dist_arr[i, j, :])

if xray:

    # definitions
    qmin = 0.1
    qmax = 8.0
    qlen = 81
    qvector = np.linspace(qmin, qmax, qlen, endpoint=True)
    elastic = False
    compton_array = x.compton_spline(atomic_numbers, qvector)

    # calculate reference IAM curve
    reference_xyz_file = "xyz/equilibrium.xyz"
    _, _, atomlist, reference_xyz = m.read_xyz(reference_xyz_file)
    reference_iam = x.iam_calc_compton(
        atomic_numbers, reference_xyz, qvector, compton_array, elastic
    )

    # read target PCD file
    # for now just create it from the xyz.
    # In general, I should read in pure data (the xyz wouldn't be known)
    target_xyz_file = "xyz/target.xyz"
    _, _, atomlist, target_xyz = m.read_xyz(target_xyz_file)
    target_iam = x.iam_calc_compton(
        atomic_numbers, target_xyz, qvector, compton_array, elastic
    )
    target_pcd = 100 * (target_iam / reference_iam - 1)
    # save to file, why not
    np.savetxt("pcd/target.dat", np.column_stack((qvector, target_pcd)))

    # calculate IAM and PCD for each frame
    iam_arr = np.zeros((qlen, nframes))
    pcd_arr = np.zeros((qlen, nframes))
    chi2_arr = np.zeros(nframes)
    inv_chi2_arr = np.zeros(nframes)
    for i in range(nframes):
        print('calculating IAM and CHI2 for frame %i' % i)
        iam_ = x.iam_calc_compton(
            atomic_numbers, xyztraj[:, :, i], qvector, compton_array, elastic
        )
        iam_arr[:, i] = iam_
        pcd_ = 100 * (iam_ / reference_iam - 1)
        pcd_arr[:, i] = pcd_
        # chi2
        chi2 = np.sum((pcd_ - target_pcd) ** 2) / qlen
        chi2_arr[i] = chi2
        inv_chi2 = 1 / chi2
        inv_chi2_arr[i] = inv_chi2
        if save_individual_pcd:
            np.savetxt(
                "pcd/%s_sample.dat" % str(i).zfill(6), np.column_stack((qvector, pcd_))
            )

# save chi2_array
np.savez("pcd/pcd_arr.npz", qvector=qvector, pcd_arr=pcd_arr)
np.savetxt("pcd/chi2_sample.dat", chi2_arr)
np.savetxt("pcd/inv_chi2_sample.dat", inv_chi2_arr)

# binning distances.
## distances i, j
for i in range(nind):
    for j in range(i + 1, nind):
        r = dist_arr[i, j, :]
        
        start = np.min(r)
        print("bin range start: %f" % start)
        end = np.max(r)
        print("bin range end: %f" % end)
        Nbins = 60
        bins = np.linspace(start, end, Nbins, endpoint=True)
        print("bin separation: %8.6f" % (bins[1] - bins[0]))
        print(bins)
        
        inds = np.digitize(r, bins)
        
        with open("r%i%i_acc.dat" % (i, j), "w") as f:
            for i in range(Nbins):
                # ignore bins with low statistics
                if len(inv_chi2[inds == i]) > 10:
                    f.write("%10.8f %10.8f \n" % (bins[i], np.max(inv_chi2[inds == i])))


print("Total time: %3.2f s" % float(default_timer() - start))
