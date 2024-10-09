# cardy

This repository contains parallelised MCMC code in C++ to sample from very high dimensional probability distributions. Our calculations can require precision much higher than the standard `long double` which is usually 80 bits. Thus this repo uses GMP and MPFR to perform arbitrary precision calculations. 

On a 64 core machine, this code is snappy enough to produce 5.4 billion single-dimensional samples per hour! [^1]

## Instructions to run on Baobab cluster at UniGE

The first step is to set the simulation parameters. Currently this is handled inside the `main.cpp` file by changing the values of `thisSimulation`. 

Then to load the required modules, run
```
module load GCC/12.3.0 GMP MPFR Boost Arrow/11.0.0 CMake
```

To build the executables, `cd` into the build directory and run (assuming `..` to the build directory contains the `CMakeLists.txt` file)
```
cmake .. && make
```

To run the simulation, use the command
```
sbatch name_of_sbatch_directive /path/to/folder
```

The `sbatch` directive file should look like this 
```
#SBATCH --job-name=mcmc
#SBATCH --partition=private-dpt-cpu
#SBATCH --time=00-03:00:00
#SBATCH --output=/path/to/data/slurm-output/slurm-outfiles/%j.out
#SBATCH --error=/path/to/data/slurm-output/slurm-errors/%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=30G

export OMP_NUM_THREADS=10;
cd /path/to/build
./cardy $@
```

The above directive is for the situation where we parallelise over 10 cores. Change as `--ntasks` and `OMP_NUM_THREADS` as desired. The code is stupidly parallelisable so to get the best performance it is advised to set `--ntasks` the same as `thisSimulation.numWalkers`.

## Footnotes
[^1]: 64 parallel random walkers run for one hour produce 54.5 million samples of a 100 dimensional distribution on a AMD EPYC-7742 2.25GHz CPU 
