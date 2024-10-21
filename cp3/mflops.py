import numpy as np

# Problem sizes and their elapsed times
problem_sizes = np.array([1024, 2048, 4096, 8192, 16384])
elapsed_times = np.array([0.00088, 0.00358, 0.01460, 0.05827, 0.23393])  # in seconds

# Calculate MFLOP/s for each problem size using the formula: MFLOP/s = (2 * n^2) / (time * 10^6)
mflops = (2 * problem_sizes ** 2) / (elapsed_times * 1e6)

print(mflops)
