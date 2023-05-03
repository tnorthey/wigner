import numpy as np
from scipy.optimize import curve_fit
import sys

# my module
import molecule

# initiate class instances
m = molecule.Molecule()
x = molecule.Xray()
sa = molecule.Simulated_Annealing()

npz_data_file = str(sys.argv[1])
# number of bins to accumulate distances into
nbins = 30
A = 100   # Strength of HO contribution


def fit_pool(npz_data_file, nbins):
    """
    reads pool.npz, then fits to target_pcd
    """

    # read pool.npz
    data = np.load(npz_data_file)  # load
    natoms = data["natoms"]
    atomarray = data["atomarray"]
    atomic_numbers = data["atomic_numbers"]
    nframes = data["nframes"]
    xyzpool = data["xyzpool"]
    distance_indices = data["distance_indices"]
    dist_arr = data["dist_arr"]
    reference_xyz = data["reference_xyz"]
    reference_iam = data["reference_iam"]
    qvector = data["qvector"]
    elastic_bool = data["elastic_bool"]
    pcd_arr = data["pcd_arr"]
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

    ### Harmonic oscillator for each distance index
    # ?
    def ho_value(A, r, r0):
        return A * (r - r0)**2
    #
    #

    # calculate CHI2 and 1/CHI2 for each structure
    chi2_arr = np.zeros(nframes)
    inv_chi2_arr = np.zeros(nframes)

    # CHD specific; notably not indices 0, 5 (the ring-opening)
    ho_indices1 = [0, 1, 2, 3, 4 ] 
    ho_indices2 = [1, 2, 3, 4, 5 ] 
    nho_indices = len(ho_indices1)

    r0_arr = np.zeros(nho_indices)
    for i in range(nho_indices):
        r0_arr[i] = np.linalg.norm(reference_xyz[ho_indices1[i], :] - reference_xyz[ho_indices2[i], :])

    total_harmonic_contrib = 0
    total_xray_contrib = 0
    print("calculating CHI2 for %i frames..." % nframes)
    for k in range(nframes):
        
        # harmonic oscillator part of chi2
        harmonic_contrib = 0
        for i in range(nho_indices):
            r = dist_arr[ho_indices1[i], ho_indices2[i], k]
            harmonic_contrib += ho_value(A, r, r0_arr[i])
        total_harmonic_contrib += harmonic_contrib

        #print('harmonic contrib: %f' % harmonic_contrib)
        chi2 = np.sum((pcd_arr[:, k] - target_pcd) ** 2) / qlen
        #print('x-ray contrib: %f' % chi2)
        total_xray_contrib += chi2
        chi2 += harmonic_contrib
        #print('chi2: %f' % chi2)
        inv_chi2 = 1 / chi2
        chi2_arr[k] = chi2
        inv_chi2_arr[k] = inv_chi2

    total_contrib = total_xray_contrib + total_harmonic_contrib
    xray_ratio = total_xray_contrib / total_contrib
    harmonic_ratio = total_harmonic_contrib / total_contrib
    print('xray contrib ratio: %f' % xray_ratio)
    print('harmonic contrib ratio: %f' % harmonic_ratio)

    ### Binning and fitting the Gaussian
    # Fitting function
    def func(xx, a, x0, sigma):
        return a * np.exp(-((xx - x0) ** 2) / (2 * sigma**2))

    # binning distances.
    ## distances i, j
    mu_arr = np.zeros((nind, nind))

    for i in range(nind):
        for j in range(i + 1, nind):
            r = dist_arr[i, j, :]
            print('distance %i %i' % (i, j))
            start = np.min(r)
            #print("bin range start: %f" % start)
            end = np.max(r)
            #print("bin range end: %f" % end)
            bins = np.linspace(start, end, nbins, endpoint=True)
            #print("bin separation: %8.6f" % (bins[1] - bins[0]))
            #print(bins)

            inds = np.digitize(r, bins)
            #print(inds)

            # put 1/chi2 into the bins
            inv_chi2_bins = np.zeros(nbins)
            for k in range(nbins):
                tmp = inv_chi2_arr[inds == k]
                # ignore bins with low statistics
                if len(tmp) > 10:
                    # in this way there can be values of 0 in the array, maybe bad for fitting?
                    # I could use append instead ?
                    ## only save the maximum 1/chi2 value in each bin
                    inv_chi2_bins[k] = np.max(tmp)

            # do I want these final dat files here? (for now ok)
            with open("results/r%i%i_acc_%i.dat" % (i, j, nframes), "w") as f:
                for k in range(nbins):
                    f.write("%10.8f %10.8f \n" % (bins[k], inv_chi2_bins[k]))

            # Executing curve_fit on noisy data
            xn = bins
            yn = inv_chi2_bins
            # popt returns the best fit values for parameters of the given model (func)
            #print('distance bins:')
            #print(xn)
            #print('1/chi2 per bin:')
            #print(yn)
            popt, pcov = curve_fit(func, xn, yn)
            # popt[0] = A,  popt[0] = mu,  popt[0] = sigma 
            print('optimised parameters (A, mu, sigma):')
            print(popt)
            mu_arr[i, j] = popt[1]

            ym = func(xn, popt[0], popt[1], popt[2])
            #print('Gaussian fit:')
            #print(ym)

    print('r0_arr')
    print(r0_arr)

    best_dist = 1e9  # start off really large
    print('Finding closest structure to fitted distances...')
    for k in range(nframes):
        dist_min = 0
        for i in range(nind):
            for j in range(i + 1, nind):
                r = dist_arr[i, j, k]
                dist_min += ( r - mu_arr[i, j] )**2  # (only carbon-carbon distances)
        dist_min_total = dist_min
        if dist_min_total < best_dist:
            best_dist = dist_min_total
            best_frame = k

    print('closest fit frame: %i' % best_frame)
    xyz_best = xyzpool[:, :, best_frame]
    m.write_xyz('best_fit.xyz', 'best', atomarray, xyz_best)
    print("Done. Wrote 'best_fit.xyz'.")


# run function
fit_pool(npz_data_file, nbins)


