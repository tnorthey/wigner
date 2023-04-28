import numpy as np
import sys

fname = str(sys.argv[1])
print("loading %s" % fname)
x = np.loadtxt(fname)
r = x[:, 0]
inv_chi2 = x[:, 1]
print(r)

start = np.min(r)
print("bin range start: %f" % start)
end = np.max(r)
print("bin range end: %f" % end)
Nbins = 60
bins = np.linspace(start, end, Nbins, endpoint=True)
print("bin separation: %8.6f" % (bins[1] - bins[0]))
print(bins)

inds = np.digitize(r, bins)

with open("out.dat", "w") as f:
    for i in range(Nbins):
        # ignore bins with low statistics
        if len(inv_chi2[inds == i]) > 10:
            f.write("%10.8f %10.8f \n" % (bins[i], np.max(inv_chi2[inds == i])))
