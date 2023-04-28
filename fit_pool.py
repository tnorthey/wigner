import numpy as np
import sys

# my module
import molecule

# initiate class instances
m = molecule.Molecule()
x = molecule.Xray()

npz_data_file = str(sys.argv[1])
# number of bins to accumulate distances into
nbins = 30

def fit_pool(npz_data_file, nbins):
    '''
    reads pool.npz, then fits to target_pcd
    '''

    # read pool.npz
    data = np.load(npz_data_file)  # load
    natoms=data["natoms"]
    atomarray=data["atomarray"]
    atomic_numbers=data["atomic_numbers"]
    nframes=data["nframes"]
    xyztraj=data["xyztraj"]
    distance_indices=data["distance_indices"]
    dist_arr=data["dist_arr"]
    reference_xyz=data["reference_xyz"]
    reference_iam=data["reference_iam"]
    qvector=data["qvector"]
    elastic_bool=data["elastic_bool"]
    pcd_arr=data["pcd_arr"]
    qlen = len(qvector)
    nind = len(distance_indices)
    
    # calculate target PCD
    ## target_pcd should be an input, but using theoretical data for now
    target_xyz_file = "target.xyz"
    _, _, atomarray, target_xyz = m.read_xyz(target_xyz_file)
    atomic_numbers = [m.periodic_table(symbol) for symbol in atomarray]
    compton_array = x.compton_spline(atomic_numbers, qvector)
    target_iam = x.iam_calc_compton(
        atomic_numbers, target_xyz, qvector, compton_array, elastic_bool
    )
    target_pcd = 100 * (target_iam / reference_iam - 1)

    # calculate CHI2 and 1/CHI2 for each structure
    chi2_arr = np.zeros(nframes)
    inv_chi2_arr = np.zeros(nframes)
    for i in range(nframes):
        print('calculating CHI2 for frame %i' % i)
        chi2 = np.sum((pcd_arr[:, i] - target_pcd) ** 2) / qlen
        inv_chi2 = 1 / chi2
        chi2_arr[i] = chi2
        inv_chi2_arr[i] = inv_chi2
    
    # binning distances.
    ## distances i, j
    for i in range(nind):
        for j in range(i + 1, nind):
            r = dist_arr[i, j, :]
            
            start = np.min(r)
            print("bin range start: %f" % start)
            end = np.max(r)
            print("bin range end: %f" % end)
            bins = np.linspace(start, end, nbins, endpoint=True)
            print("bin separation: %8.6f" % (bins[1] - bins[0]))
            print(bins)
            
            inds = np.digitize(r, bins)
            print(inds)
            
            # do I want these final dat files here? (for now ok)
            with open("r%i%i_acc_%i.dat" % (i, j, nframes), "w") as f:
                for k in range(nbins):
                    # ignore bins with low statistics
                    tmp = inv_chi2_arr[inds == k]
                    if len(tmp) > 10:
                        f.write("%10.8f %10.8f \n" % (bins[k], np.max(tmp)))


# run function
fit_pool(npz_data_file, nbins)

