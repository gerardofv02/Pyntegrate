import numpy as np

arr1 = np.load("input_array_bam_after.npz.npy")
arr2 = np.load("input_array_bw_after.npz.npy")

print(np.array_equal(arr1, arr2))

##all cases false
