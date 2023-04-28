import numpy as np
import sys

# my module
import molecule

# initiate class instances
m = molecule.Molecule()
x = molecule.Xray()

# calculate target PCD file
target_pcd = 100 * (target_iam / reference_iam - 1)

def fit_pool(target_pcd, nbins):
    '''
    reads pool.npz, then fits to target_pcd
    '''

    # read pool.npz
    data = np.load("pool.npz")
    natoms = data["natoms"]
    atomarray = data["atomarray"]
    atomic_numbers = data["atomic_numbers"]
    nframes = data["nframes"]
    xyztraj = data["xyztraj"]
    pcd = ### NOT FINISHED!!! ####
    
    # calculate CHI2 and 1/CHI2 for each structure
    chi2_arr = np.zeros(nframes)
    inv_chi2_arr = np.zeros(nframes)
    for i in range(nframes):
        print('calculating CHI2 for frame %i' % i)
        chi2 = np.sum((pcd_ - target_pcd) ** 2) / qlen
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
            
            # do I want these final dat files here? (for now ok)
            with open("r%i%i_acc.dat" % (i, j), "w") as f:
                for i in range(nbins):
                    # ignore bins with low statistics
                    if len(inv_chi2[inds == i]) > 10:
                        f.write("%10.8f %10.8f \n" % (bins[i], np.max(inv_chi2[inds == i])))


