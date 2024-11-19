
```markdown
# Sobel Filter GPU Performance Study

This project explores GPU programming techniques for image processing by implementing a Sobel filter for edge detection using three different parallel computing models: OpenMP on the CPU, CUDA on the GPU, and OpenMP with device offload on the GPU. The following instructions guide you through setting up, building, and running the code on the Perlmutter system.

## Prerequisites

To prepare the environment for GPU execution, ensure the following modules and environment variables are configured:

```bash
module load PrgEnv-nvidia 
export CC=cc
export CXX=CC
```

## Build Instructions

1. **Clone the repository and create a build directory:**
    ```bash
    cd $SCRATCH/your-code-folder
    mkdir build
    cd build
    cmake ../
    make
    ```

## Running the Code

The codebase includes three main implementations:
1. `sobel_cpu.cpp` – CPU implementation with OpenMP
2. `sobel_gpu.cu` – GPU implementation using CUDA
3. `sobel_cpu_omp_offload.cpp` – GPU implementation with OpenMP offload

Navigate to your `$SCRATCH` directory before running any of the commands below.

### 1. Running the CUDA Code

To run the CUDA code (`sobel_gpu.cu`) with different configurations of threads per block and block counts, use the following SLURM script to automatically iterate over configurations and gather profiling metrics.

Create a file named `sobel_gpu_run.sh` with the following content:

```bash
#!/bin/bash 
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -G 1
#SBATCH -q regular
#SBATCH -t 01:00:00
#SBATCH -A m3930
#SBATCH -J gpu-job
#SBATCH --output=gpu-job.o%j
#SBATCH --error=gpu-job.e%j

# Change to the build directory
cd $SCRATCH/build

# Define the configurations
threads_per_block=(32 64 128 256 512 1024)
num_blocks=(1 4 16 64 256 1024 4096)

# Run the code with different configurations and collect performance metrics
for tpb in "${threads_per_block[@]}"
do
    for nb in "${num_blocks[@]}"
    do
        echo "Running with $tpb threads per block and $nb blocks"
        ncu --set default --section SpeedOfLight --section Occupancy --section LaunchStats --section SourceCounters --metrics smsp__cycles_active.avg.pct_of_peak_sustained_elapsed,dram__throughput.avg.pct_of_peak_sustained_elapsed,gpu__time_duration.avg ./sobel_gpu $nb $tpb
    done
done
```

Submit the job using:
```bash
sbatch sobel_gpu_run.sh
```

### 2. Running the CPU Code

For the CPU implementation, specify the concurrency level directly as an argument to `sobel_cpu`:

```bash
./sobel_cpu <concurrency_level>
```

For example, to run with a concurrency level of 16 threads:
```bash
./sobel_cpu 16
```

This command will execute the CPU-based Sobel filter with the specified number of threads.

### 3. Running the OpenMP Offload Code

Since the OpenMP offload code is designed for a single configuration, you can use the following command directly to run the offload version and collect performance metrics:

```bash
cd $SCRATCH/build
ncu --set default --section SpeedOfLight --section Occupancy --section LaunchStats --section SourceCounters --metrics smsp__cycles_active.avg.pct_of_peak_sustained_elapsed,dram__throughput.avg.pct_of_peak_sustained_elapsed,gpu__time_duration.avg ./sobel_cpu_omp_offload
```

## Output

- Each run collects data on:
  - **Runtime** – Measured in milliseconds.
  - **Achieved occupancy** – Reported as a percentage.
  - **Percentage of peak sustained memory bandwidth** – Reported as a percentage.
