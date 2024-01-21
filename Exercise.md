# Running molecular dynamics simulations with GROMACS on LUMI

## Introduction
==CB:Write an introductiory paragraph about the aim of the tutorial and some information about LUMI?==

==CB: I feel like goals shouldn't be formulated in present continuous form so I changed that. Feel free to change back if you don't agree! AV Usually ILOs are in future tense, goals/purpuses in to + infinitive==
:::success 
:dart: **Learning goals**
* Get familiar with the GROMACS tools used in the exercises.
* Understand common features of the `mdrun` command line.
* Understand key parts of the `mdrun` log file structure.
:::

Authors: Szilárd Páll and Andrey Alekseenko (adaped from [doi:10.6084/m9.figshare.22303477](https://doi.org/10.6084/m9.figshare.22303477))

Software: GROMACS-2023.3

[//]: # (### Performance in MD)
[//]: # (Molecular dynamics simulations do the same calculation repeatedly, typically for a large number of time-steps. )

### The GROMACS simulation engine 

GROMACS is a molecular simulation package which comes with a molecular dynamics simulation engine (`mdrun`), a set of analysis tools, and the gmxlib Python API. 

GROMACS is highly flexible and can be built in various ways depending on the target harware architecture and the parallelization features enabled. The GROMACS features and dependencies enabled at compile-time are shown in the _version header_ which is listed by invoking the `gmx -version` command as well as at the top of the simulation log outputs. All functionalities of the GROMACS package, the simulation engine and tools, are provided in the `gmx` program through subcommands. The program can have suffixes, e.g. MPI builds are typically installed as `gmx_mpi`.

==CB: Maybe add a link to the GROMACS manual or installation instructions here?==
==CB: Maybe write something like: "In this tutorial, we will use a version of GROMACS that has already been built on the LUMI-G cluster, but if you wish to install GROMACS on your own system, instructions for many different hardware configurations are available at [GROMACS documention]."==

#### GROMACS parallelization overview

Parallelization of MD simulation requires _expressing_ concurrent work (multiple computations happening at the same time) and _exposing_ it using an implementation with the help of a parallel programming model. To express concurrency within a single simulation in GROMACS we can divide the work using data (e.g. spatial decomposition algorithms), task (e.g. rank specialization for the "separate PME ranks" feature), or ensemble decomposition. The exposed concurrent work can then be mapped to various processing units, like CPU cores or GPU accelerators. GROMACS realies on a hierarchical heterogeneous parallelization using MPI, OpenMP multi-threading, CUDA/SYCL/OpenCL for asynchronous GPU execution, and SIMD for low-level CPU and GPU algorithms.

The data parallelism is used for implementing spatial decomposition (that consist of dividing the production system into subsystems that are as independent as possible) and takes place
across MPI ranks using multi-threading on CPUs and fine-grained SIMD-style algorithms (Single Instruction Multiple Data). At the same time, task parallelism is at the heart of the heterogeneous GPU engine and it is also what enables scaling the PME algorithm efficiently by employing [rank specialization](https://manual.gromacs.org/documentation/current/user-guide/mdrun-performance.html#separate-pme-ranks).

MD simulation studies can be classified into two major flavors: a single (or a few) long trajectory or larger set of trajectories. Due to the timescale of some biological processes,
a single/few very long trajectories might not be enough (and inefficient) to sample the conformational space. Then, an alternative is to use an ensemble of simulations. 

A wide range of algorithms, from free energy perturbation to replica exchange to the accelerated weight histogram method (AWH), rely on (or require) multiple MD simulations which form an _ensemble_. An _ensemble simulation_ refers to a set of simulations, where each individual simulation is referred to as _ensemble member_ (called "replica" in replica-exchange and walker in AWH). These algorithms provide a source of concurrent work, which simulation workflows can use to parallelize over, and require different levels of coupling between the ensemble members. E.g., standard free-energy calculations (with a pre-determined simulation length) requires no communication across the ensemble members, whereas replica-exchange and AWH require exchange of information at regular time intervals. The latter class of methods is referred to as _coupled_ ensembles. Depending on the frequency of data exchange, ensembles can be _weakly_ or _strongly_ coupled (with infrequent or frequent data exchange, resp.). Coupled ensembles are more performance sensitive, hence more prone to be influenced by imbalance (e.g. member simulations running with different throughput). The stronger the coupling the more sensitive the ensemble simulation is to performance bottlenecks.


#### The `mdrun` simulation tool

In GROMACS, the primary way to run simulations across CPUs and GPUs is to use the command line program `mdrun`. The simulation tool `mdrun` can be invoked as a subcommand of the main program, e.g. `gmx mdrun`. The `mdrun` functionalities available in a specific build depend on the options GROMACS was configured with and can be seen in the version header.

The following list contains the central performance-related command line options used in this tutorial:
* `-g LOGFILE` set a custom name for the log file (default `md.log`);
* `-pin on` enable `mdrun` internal thread affinity setting (might override externally set affinities). Note on LUMI externally set affinities are recommended or [the LUMI documentation](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/distribution-binding/);
* `-tunepme`/`-notunepme` enable PME task load balancing;
* `-nsteps N` set the number of simulations steps for the current run to `N` (N=-1 means infinitely long runs, intended to be combined with `-maxh`);
* `-maxh H` stop the simulation after `0.99*H` hours;
* `-resethway` reset performance counters halfway through the run;
- `-nb`/`-pme`/`-bonded`/`-update` task assignment options used to select tasks to run on either CPU or GPU. 
- `-npme N` set the number of separate ranks to be used for PME (N=-1 is a guess)

==AV to AA/SP what about to add the following sentences:  "Note that some performance features can be triggered using environment variables. See examples in [exercise 3.2](https://hackmd.io/qvLmXFLCQGScOHhdjS5uQw?both#Exercise-32-Separate-PME-rank-with-direct-GPU-communication)"" == 

For further information on the `mdrun` simulation tool command line options and features, see the [online documentation](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html).


#### The `mdrun` log file

==CB: The figures need more explaining. There are also multiple versions of the same figures.==

The log file of the `mdrun` simulation engine contains extensive information about the GROMACS build, hardware detected at runtime, complete set of simulation setting, diagnostic output related to the run and its performance, as well as physics and performance statistics.

The version header in a GROMACS `mdrun` log ==AV to AV change picture to 1MPI/1GPU==:
![Picture1bisbis](https://hackmd.io/_uploads/rycDUsYdT.jpg)

The hardware detection section of the `mdrun` log:
![Picture2](https://hackmd.io/_uploads/SkM9l9cua.png)

![](https://i.imgur.com/68AWSW1.png)

The performance accounting in a multi-rank simulation:

![Picture3](https://hackmd.io/_uploads/S1BAS95_a.jpg)![Picture4](https://hackmd.io/_uploads/BkMr_cc_p.png)



![](https://i.imgur.com/PUlYmaR.png)


## Simulation input system

In the following exercises, we will use two different simulation systems: 
* large-sized [satellite tobacco mosaic virus, STMV](https://en.wikipedia.org/wiki/Tobacco_virtovirus_1) (~ 1 milion atoms) system solvated in a box of TIP3P water molecules, using the CHARMM27 force field.

[//]: # (Comment 1066628 atoms, 2 fs time-step, 1.2 nm cut-offs, h-bond constraints, 0.15 nm PME grid spacing, NVT ensemble.)

* medium-sized [aquaporin membrane protein](https://en.wikipedia.org/wiki/Aquaporin), a tetrameric ion channel (~110000 atoms) embedded in a lipid bilayer and solvated in a box of TIP3P water using the CHARMM36 force field. We will use the Accelerated Weight Histogram (AWH) algorithm with 32 walkers to enhance sampling of the simulations
==CB:What specifically will be enhanced in AWH?==

[//]: # (Comment 2.5 fs time-step, 1.2 nm cut-offs, 0.1125   nm PME grid spacing, h-bond constraints, NPT ensemble.)

Both systems have previously been used to benchmark GROMACS heterogenous parallelization and acceleration (https://doi.org/10.1063/5.0018516). 

The simulation input files (tpr) can be obtained from:

==AV comment to SP: not sure that they are consistent with LUMI files [STMV](https://kth-my.sharepoint.com/:u:/g/personal/pszilard_ug_kth_se/EZUGKtf30vlOtMtjyjpjGVwBGOXCNEIkcp2BFoZDzme5gQ) [Aquaporin](https://kth-my.sharepoint.com/:u:/g/personal/pszilard_ug_kth_se/Eawzt6R9hrlPjub4ZXBYOVMBPKdx0ZNRn0G4GSLG6ndYiA)==

* [workshop input files](https://a3s.fi/gmx-lumi/workshop-input-files.tar.gz)
* or in the folder `/projappl/project_465000934` on LUMI

## 1. Running your first jobs on LUMI-G

In this first exercise, we will submit our initial jobs on LUMI-G and explore key features and peculiarities of the LUMI system, scheduler, and GROMACS' `mdrun` simulation tool. As simulation system we use the satellite tobacco mosaic virus, STMV. You find the input file in the folder `/projappl/project_465000934/stmv` on LUMI 
==AV to check that the directory path is correct)==

We will start with a basic job submission script (batch script) and successively build on it to explore how to correctly request resources using the SLURM job scheduler and to finally arrive to have a script that correctly requests resources on LUMI-G nodes.

:::info
:bulb: The LUMI-G hardware partition consists of 2978 nodes with 4 AMD MI250x GPUs and a single 64 cores AMD EPYC "Trento" CPU. Each MI250x GPU is a multi-chip module with two GPU dies named by AMD Graphics Compute Die (GCD). For furtherdetails on the LUMI architecture, see [the LUMI documentation](https://docs.lumi-supercomputer.eu/hardware/lumig/) or [how to use GROMACS on LUMI](https://docs.csc.fi/apps/gromacs/).
:::

:::success
:dart: **Learning goals**
* Know how to submit GROMACS jobs.
* Get familiar with common SLURM scheduling and `mdrun` command line options.
* Get familiar with the `mdrun` console and log outputs.
* _Bonus_: Understand the impact of using multiple CPU cores/threads on parts of the MD computation.
:::

:::info
:bulb: **The GROMACS log file** Before starting, take a look at the [introduction on the GROMACS`mdrun` simulation tool](https://hackmd.io/qvLmXFLCQGScOHhdjS5uQw#The-mdrun-simulation-tool) and at the description of log files.
:::

### Exercise 1.1: Launching a first GROMACS simulation on LUMI-G
Now, we will launch a first test simulation on LUMI-G. Make sure you have copied the necessary topol.tpr file into a working directory of your choice (this is also where your output files will end up, so it's good to keep it organized!), and then create a batch file (with suffix .sh) with the following content:

```bash=
#!/bin/bash
#SBATCH --partition=small-g            # partition to use
#SBATCH --account=project_465000934    # project for billing
#SBATCH --reservation=gromacs_wednesday # workshop reservation 
#SBATCH --time=00:10:00                # maximum execution time of 10 minutes
#SBATCH --nodes=1                      # we use 1 node
#SBATCH --ntasks-per-node=1            # 1 MPI rank
#SBATCH --cpus-per-task=1              # number cpus-per-task

module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun gmx_mpi mdrun -g ex1.1_${SLURM_NTASKS}x${OMP_NUM_THREADS}_jID${SLURM_JOB_ID} -nsteps -1 -maxh 0.017 -resethway -notunepme

```
Note the four benchmarking flags used above (`-nsteps`, `-maxh`, `-resethway` and `-[no]tunepome`) [are described in the introduction](https://hackmd.io/qvLmXFLCQGScOHhdjS5uQw#The-mdrun-simulation-tool).
We use time limited runs here by passing `-maxh 0.017` which sets the run time limit in hours (~one minute); we do that as the simulation throughput significantly changes and we want to avoid either very long wait time or unreliable benchmark measurements due to just a few seconds of runtime.

==CB: Here I would add something like: "To submit the job, use the `sbatch` command and wait until it finishes (it should not take longer than 10 minutes!)."==

:::warning
* Take a look at log file (.log) and find the hardware detection and performance table. What are the resources detected and used? 
:::

Now try to enable multithreading. To do that, we need to request multiple CPU cores. Edit the job script, change the number of CPU cores and submit a new job.

```bash=
#!/bin/bash
#SBATCH --partition=small-g            # partition to use
#SBATCH --account=project_465000934    # project for billing
#SBATCH --reservation=gromacs_wednesday # reservation for 24 January 
#SBATCH --time=00:10:00                # maximum execution time of 10 minutes
#SBATCH --nodes=1                      # we run on 1 node
#SBATCH --ntasks-per-node=1            # 1 MPI rank
#SBATCH --cpus-per-task=...            # number cpus-per-task

module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun gmx_mpi mdrun -g ex1.1_${SLURM_NTASKS}x${OMP_NUM_THREADS}_jID${SLURM_JOB_ID} -nsteps -1 -maxh 0.017 -resethway -notunepme

```

:::warning
* Compare the log file hardware detection and performance table of the two above runs. What has changed in terms of resources detected and used? Is there anything still missing? 
:::

LUMI-G has relatively few CPUs cores per GPU, so making the best use of these is important and can have a strong impact on performance. We will explore this further in Exercise 1.3.

### Exercise 1.2: Launching a simple GROMACS GPU run

==CB: I don't understand how it can be seen in the previous exercise that we need to request a GPU. We just want a GPU? AV I edit below based to CB comment AA/SP check that it is what you want==
As note in the previous exercise, GPU are not detected in the log file, but a LUMI-G node has 8 graphical cards. Now we learn how to request GPUs. Use the job script below to submit a job using one GCD (in particular, this is usefull for jobs which do not use the full node).

== AV to SP and AA do we use both GCD and GPU?== 

```bash=
#!/bin/bash
#SBATCH --partition=small-g            # partition to use
#SBATCH --account=project_465000934    # project for billing
#SBATCH --reservation=gromacs_wednesday # reservation for 24 January 
#SBATCH --time=00:10:00                # maximum execution time of 10 minutes
#SBATCH --nodes=1                      # we run on 1 node
#SBATCH --ntasks-per-node=1            # 1 MPI rank
#SBATCH --cpus-per-task=7              # number cpus-per-task
#SBATCH --gpus-per-task=1              # New line! Get 1 GPU device  

module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun gmx_mpi mdrun -g ex1.2_${SLURM_NTASKS}x${OMP_NUM_THREADS}_jID${SLURM_JOB_ID}
-nsteps -1 -maxh 0.017 -resethway -notunepme

```

:::warning
* Look at the log file hardware detection and performance table: what are the resources detected and used?
* Compera the log file of Ex 1.1 and 1.2. Has the performance changed?
:::

### _Bonus_ Exercise 1.3: Explore the use of CPUs and OpenMP multi-threading

==AA to SP/AA I am confuse, in the example script below you ask gpus, but you wrote that we do not use GPU.. I guess the script is wrong. Or?==

In this exercise, we will use only the CPUs of the LUMI-G nodes to explore how the different computational tasks perform with OpenMP multi-threading. 

```bash=
#!/bin/bash
#SBATCH --time=00:10:00        # maximum execution time of 10 minutes
#SBATCH --partition=dev-g
#SBATCH --account=project_...  # Project for billing
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-task=1      # New line! Get a GPU

module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun gmx_mpi mdrun -nsteps -1 -maxh 0.017 -resethway -notunepme \
    ex1.3_1x${SLURM_CPUS_PER_TASK}_jID${SLURM_JOB_ID}
```
==CB: Here, instructions come after the script. Previously, they came before. Maybe good to be consistent?==
* Modify the script varying the number of CPU cores used (`--cpus-per-task`) and submit runs with each new configuration.
* Look at the `mdrun` log file output (the files will be named `ex1.3_1xN_jIDXXXXXX.log`).

:::warning
* How does the absolute performance (ns/day) change when increasing the number of cores used?
* How does the wall-time of various computations changes with the thread count (e.g. "Force", "PME mesh", "Update" tasks).

[//]: # (Comment: * _Bonus_: )
:::

:::info
:bulb: `SBATCH` arguments provided in the job script header can also be passed on the command line (e.g. `sbatch --cpus-per-task N`) overriding the setting in the job script header. Doing so can allow varying submission parameters without having to edit the job script.
:::

:::spoiler Help with the solution
Sample log files for the exercise session.: 
* [ex 1.1](https://github.com/Lumi-supercomputer/gromacs-on-lumi-workshop/tree/main/Exercise-1.1/STMV)
* [ex 1.2](https://github.com/Lumi-supercomputer/gromacs-on-lumi-workshop/tree/main/Exercise-1.2/STMV)
:::

## 2. GPU accelerated simulations
:::success
:dart: **Learning goals**
* Understand how the GROMACS heterogeneous parallelization allows moving tasks between CPU and GPU and how that impacts performance.
* Understand the difference between force-offload and GPU-resident modes.
* _Advanced_: Explore the effects of load balancing.
:::

The GROMACS MD engine uses heterogeneous parallelization which can flexibly utilize both CPU and GPU resources. As discussed in the lecture, there are two offload modes:
* In the _force offload_ mode, some or all forces are computed on the GPU, but are transferred to the CPU every iteration for integration;
* In the GPU-resident mode, the integration happens on the GPU allowing the simulation state to reside on the GPU for tens or hundreds of iterations. 
Further details can be found in the [GROMACS users guide](https://manual.gromacs.org/current/user-guide/mdrun-performance.html#running-mdrun-with-gpus) and [DOI:10.1063/5.0018516](https://aip.scitation.org/doi/full/10.1063/5.0018516).


In the following exercises, we will learn how moving tasks between the CPU and GPU impacts performance. As simulation system we use the satellite tobacco mosaic virus, STMV. You can take topol.tpr file from the exercise-1 serie (or copy it from the folder `/projappl/project_465000934/stmv` on LUMI) 
==AV to check that the directory path is correct==


We will be using LUMI-G GPU nodes for submitting single-GPU device jobs (hence using one of the eight in the full compute node); for further details on the architecture and usage see the [GPU nodes - LUMI-G](https://docs.lumi-supercomputer.eu/hardware/lumig/).


### Exercise 2.1: GPU offloading force computations

The tasks corresponding to the computation of bonded, short and long-range non-bonded forces can be offloaded to a GPU in GROMACS. The assignment of these tasks is controlled by the following `mdrun` command line options:
* (short-range) nonbonded: `-nb ASSIGNMENT`
* particle mesh Ewald: `-pme ASSIGNMENT`
* bonded: `-bonded ASSIGNMENT`

The possible "`ASSIGNMENT`" values are `cpu`, `gpu`, or `auto`.

We use one GPU with CPU cores (an eighth of the node) in a simulation and assess how the performance changes with offloading different force calculations. As a baseline, launch a run first with assigning all tasks to the CPU (as below).

```bash=
#!/bin/bash
#SBATCH --partition=small-g            # partition to use
#SBATCH --account=project_465000934    # project for billing
#SBATCH --reservation=gromacs_thursday # reservation for 25 January
#SBATCH --time=00:10:00                # maximum execution time of 10 minutes
#SBATCH --nodes=1                      # we run on 1 node
#SBATCH --gpus-per-task=1              # we use 1 GPU device 
#SBATCH --ntasks-per-node=1            # 1 MPI rank
#SBATCH --cpus-per-task=7              # number cpus-per-task

module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun gmx_mpi mdrun -nb gpu -pme gpu -bonded gpu -update cpu \
   -g ex2.1_${SLURM_NTASKS}x${OMP_NUM_THREADS}_jID${SLURM_JOB_ID} -nsteps -1 -maxh 0.017 -resethway -notunepme

```
* Next submit jobs by incrementally offloading various force tasks ( non-bonded `-nb`, PME `-pme`, bonded `-bonded`) to the GPU.

:::warning 
* How does the performance (ns/day) change with offloading more tasks?
* Look at the performance table in the log and observe how the fraction wall-time spent in the tasks left on the CPU change. 

**Note** that the log file performance report will only contain timings of tasks executed on the CPU, not those offloaded to the GPU, as well as timings of the CPU time spent launching GPU work as well as waiting for GPU results.
:::
:::info
:bulb: The mdrun option `-maxh 0.017` sets the run time limit in hours (~one minute); we do that as the simulation throughput significantly changes and we want to avoid either very long wait time or unreliable benchmark measurements due to just a few seconds of runtime.
:::

### Exercise 2.2: GPU-resident mode

Continuing from the previous exercise, we will now explore using the GPU-resident mode.

```bash=
#!/bin/bash
#SBATCH --partition=small-g            # partition to use
#SBATCH --account=project_465000934    # project for billing
#SBATCH --reservation=gromacs_thursday # reservation for 25 January
#SBATCH --time=00:10:00                # maximum execution time of 10 minutes
#SBATCH --nodes=1                      # we run on 1 node
#SBATCH --gpus-per-task=1              # we use 1 GPU device 
#SBATCH --ntasks-per-node=1            # 1 MPI rank
#SBATCH --cpus-per-task=7              # number cpus-per-task

module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun gmx_mpi mdrun -nb gpu -pme gpu -bonded gpu -update gpu \
     -g ex2.2_${SLURM_NTASKS}x${OMP_NUM_THREADS}_jID${SLURM_JOB_ID} -nsteps -1 -maxh 0.017 -resethway -notunepme
     
```

* Submit a fully offloaded GPU-resident job using the `-update gpu` option (as above).
* Since we moved most computation to the GPU, the CPU cores are left unused. The GROMACS heterogeneous engine allows moving work back to the CPU. Now, let's try to utilize CPU cores for potential performance benefits. First, try moving the PME tasks back to the CPU, then the bonded tasks.

[\\]: # (Next let's try adding more cores to the run by specifying `--cpus-per-task=` and `OMP_NUM_THREADS=` in the job script.)

:::warning 
* How does the GPU-resident mode perform compared to the best performing force-offload run from ex 2.1?
* How did the performance change when force tasks were moved back to the CPU?
* _Bonus_: Enable (PME) task load balancing by replacing `-notunepme` with `-tunepme`. Assign the PME task to the CPU and GPU and observe how the performance changes compared to the earlier similar run without load balancing.
:::

==TODO: add -nstlist tuning exercise / question== 

:::spoiler help with the results
Sample log files for the exercise session.: 
* [ex 2.1](https://github.com/Lumi-supercomputer/gromacs-on-lumi-workshop/tree/main/Exercise-2.1/STMV)
* [ex 2.2](https://github.com/Lumi-supercomputer/gromacs-on-lumi-workshop/tree/main/Exercise-2.2/STMV)
:::

## 3. Scaling GROMACS across multiple GPUs


:::success
:dart: **Learning goals**
* Understand how task- and domain decomposition is used in the `mdrun` simulation engine.
* Explore the multi-GPU strong scaling of GROMACS simulations and the effect of using different decomposition and communication schemes.

:::

In these exercises, we will explore the use of multiple GPUs (on one or several GPU-enabled nodes) to improve the simulation performance. As simulation system we use the satellite tobacco mosaic virus, STMV. You can take topol.tpr file from the previuos exercise serie (or copy it from the folder `/projappl/project_465000934/stmv` on LUMI) 
==AV to check that the directory path is correct==

### Exercise 3.1: Separate PME rank

In the previuos exercise serie, we have learnt how to offload some _tasks_ from CPU to GPU: short-range nonbonded, PME and bonded forces, update and constraints. When scaling across **two** GPU devices, the same approach can be applied. 

==CB: Nice with these explanations!==
:::info
Recall from the lecture the notion of _MPI rank_: MPI ranks coordinate their work by exchanging data over a special protocol and, unlike OpenMP threads, can run on different _nodes_ in the cluster. Each MPI rank can only use one GPU, but different MPI ranks can use different GPUs.
:::

With two MPI ranks, `mdrun` can do non-bonded, bonded, and integration tasks on the first rank and PME on the second rank. This is a basic version of run mode called _separate PME rank_. Furthermore, each task can be run either on CPU or GPU; when running on CPU, the work can be distributed across multiple OpenMP threads within a rank.
==CB: "...when running on CPU, the work can be distributed across multiple OpenMP threads within a rank." What about the GPU?==

:::info
Recall from the lecture on LUMI architecture that there is an intricate interconnection between CPUs and GPUs. To make sure the code runs optimally, always use ==TODO== wrapper script, which pins CPU and GPU tasks to the devices in the optimal fashion for LUMI.
:::

Try running on two GPUs using two ranks, one particle-particle (PP) and one PME:

==AV to AA Note that the batch scripts used for 3.1 and 3.2,  do not include `--exclusive`.  I have copied the window in 3.3 where we have --exclusive. ==

:::warning 
⚠️ ==TODO== The script below uses `--exclusive` flag, reserving the whole node, yet uses only two GPUs out of eight. This makes performance measurements from short runs more predictable, but should not be used for long runs as it wastes resources. 
:::


```bash=
#!/bin/bash
#SBATCH --partition=small-g            # partition to use
#SBATCH --account=project_465000934    # project for billing
#SBATCH --reservation=gromacs_thursday # reservation for 25 January
#SBATCH --time=00:10:00                # maximum execution time of 10 minutes
#SBATCH --nodes=1                      # we run on 1 node
#SBATCH --gpus-per-task=..             # fill in number of  GPU devices 
#SBATCH --ntasks-per-node=...          # fill in number of MPI rank
#SBATCH --cpus-per-task=7              # number cpus-per-task

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun gmx_mpi mdrun -npme 1  \
                 -nb gpu -pme gpu -bonded gpu -update gpu \
                 -g ex3.1_${SLURM_NTASKS}x${OMP_NUM_THREADS}_jID${SLURM_JOB_ID} \
                 -nsteps -1 -maxh 0.017 -resethway -notunepme

```

==TODO:Check that non-full-node runs behave consistently==

* Try to changing which tasks are offloaded to GPUs and CPUs by varying `-bonded` and `-update` flags (keeping `-nb gpu -pme gpu`). 

:::warning 
* Look at the absolute performance in the log files.
* Compare the log files when different tasks are offloaded to GPUs and CPUs. Is the GPU-resident or force-offload mode more efficient?
:::

### Exercise 3.2: Separate PME rank with direct GPU communication.

[\\]: # (In the previous exercise, you might have noticed that a lot of time is spent in communication.)

By default, GROMACS uses _staged communication_: data is copied from GPU to CPU on one rank, then sent between CPUs, and finally copied to the target GPU. On LUMI, it is more efficient to use _direct GPU communication_ (also called "GPU-aware"). This is achieved by setting the following environment variables:

```bash
export MPICH_GPU_SUPPORT_ENABLED=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1
```
* Add these three variables to the script from Exercise 3.1 (before calling `srun`) and repeat the runs of the previous exercise, but now with GPU-direct communication enabled. 

:::warning 
* Look at how the absolute performance change in the log files
* Compare the performance counters in log files in this and previous exercise. What do you observed? ==TODO==: Add hints.
* _Bonus:_ Is the GPU-resident or force-offload mode more efficient now?
:::


### Exercise 3.3: Domain-decomposition with a separate PME rank


Distributing tasks between GPUs only allows limited parallelism. If we wish to scale to more than two GPUs, _domain decomposition_ should be employed. 

==CB: I found some sentences below a bit unclear. Suggestion for reformulation: "Most interactions in molecular dynamics are short-range which generally don't interfere with the domain decomposition approach. However, long-range electrostatic interactions are less amenable to decomposition and therefore require some special attention."==

:::info
Recall from the lecture that the domain decomposition allows _spatially_ decomposing simulation data into _domains_. Most interactions in molecular dynamics are short-range except long-range electrostatics. This is less amenable to decomposition. Typically, when scaling over _N_ ranks (1 GPU per rank, _N_ > 2), we use one rank to compute long-range electrostatics (PME) for the whole system, and use domain decomposition to distribute the other tasks between the remaining _N_-1 ranks.
:::

Now we try running on 4, 6, 8 ranks (including one separate PME rank) by changing the values in the script below. To  dedicate one rank for PME we use the option `-npme 1`.
==AV to AA: If you want gpus-per-node=ntask-per-node, good specify here again (e.i something like that "ensure that you ask a number of GPUs equal to MPI ranks"== 


```bash=
#!/bin/bash
#SBATCH --partition=small-g            # partition to use
#SBATCH --account=project_465000934    # project for billing
#SBATCH --reservation=gromacs_thursday # reservation for 25 January
#SBATCH --exclusive                    # new! to reserve the whole node
#SBATCH --time=00:10:00                # maximum execution time of 10 minutes
#SBATCH --nodes=1                      # we run on 1 node
#SBATCH --gpus-per-node=..             # fill in number of  GPU devices per node
#SBATCH --ntasks-per-node=...          # fill in number of MPI rank per node (max 8)

module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu
source ${GMXBIN}/lumi-affinity.sh

export OMP_NUM_THREADS=7

export MPICH_GPU_SUPPORT_ENABLED=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1

srun --cpu-bind=${CPU_BIND} ./select_gpu \
     gmx_mpi mdrun -npme 1 -nb gpu -pme gpu -bonded gpu -update gpu \
                   -g ex3.3_${SLURM_NTASKS}x${OMP_NUM_THREADS}_jID${SLURM_JOB_ID} \
                   -nsteps -1 -maxh 0.017 -resethway -notunepme
```
==AV to AA: I suggest to add this sentence here "The script below uses `--exclusive` flag, to reserve the whole node. This makes performance measurements from short runs more predictable, but should not be used for long runs as it wastes resources."== 

==AV to AA: to be on the save side, please check. Below your old example script. Above example script from CSC. Please delete the script below if everything is cover by the above one.==

```bash=
#!/bin/bash     NOT LUMI SCRIPT
#
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=...      # Number of ranks, max. 8 per node
#SBATCH --cpus-per-task=7
#SBATCH --gres=gpu:...    # Should match the number of tasks

module load TODO
export OMP_NUM_THREADS=5
export MPICH_GPU_SUPPORT_ENABLED=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1

srun launcher_script_todo.sh gmx_mpi mdrun \
    -npme 1 \ # We always dedicate one rank for PME
    -nb gpu -pme gpu -bonded cpu -update cpu \
    -g ex3.3_${SLURM_NTASKS}x${OMP_NUM_THREADS}_jID${SLURM_JOB_ID} \
    -nsteps -1 -maxh 0.017 -resethway -notunepme
```


:::warning 
* Look at the absolute performance in the log file.
* How does absolute performance scale with increasing number of GPUs?
* Bonus: Look at the values of nstlist in the log files. What do you observe?  
:::

:::spoiler Help with the solution
Sample log files for the exercise session.: 
* [ex 3.1](https://github.com/Lumi-supercomputer/gromacs-on-lumi-workshop/tree/main/Exercise-3.1/STMV)
* [ex 3.2](https://github.com/Lumi-supercomputer/gromacs-on-lumi-workshop/tree/main/Exercise-3.2/STMV)
* [ex 3.3](https://github.com/Lumi-supercomputer/gromacs-on-lumi-workshop/tree/main/Exercise-3.3/STMV)
:::

### Exercise 3.4: Domain-decomposition with a PME decomposition (advanced)

==TODO== Decide whether a separate module is needed and whether scaling works at all

In the previous exercises, we explored scaling of GROMACS simulation across up to 8 GPUs using task-decomposition to assign long-range electrostatic work to one GPU, and domain-decomposition to distribute the short-range work between the remaining GPUs. Past a certain number of ranks, PME work running on a single device might become the bottleneck. In this case, _PME decomposition_ must be used. To use it, you can either:
- set `-npme` to a value greater than 1, in which case the specified number of ranks will be dedicated to PME work, and the rest to the short-range work, or
- set `-npme 0`, in which case all ranks will do short-range and long-range work.

==CB:"all ranks will do **both(?)** short-range and long-range work"==

Try running the simulation on two nodes (16 GPUs), with 1, 2, and 4 dedicated PME ranks, as well as with PME decomposition over all ranks (`-npme 0`):

==AV to AA we do not have example script, neither test in  CSC repo. Shall we tell the attendees to modify one of the script above (if we are sure that it works on 2 nodes).==

```bash=
#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=7
#SBATCH --gres=gpu:16

module load TODO
export OMP_NUM_THREADS=5
export MPICH_GPU_SUPPORT_ENABLED=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1
export GMX_GPU_PME_DECOMPOSITION=1 # Needed to enable PME decomposition

NPME=... # Set the number of dedicated PME ranks (or 0 to run PME on all PP ranks)

srun launcher_script_todo.sh gmx_mpi mdrun \
	-npme ${NPME} \ # Decompose PME over two ranks
	-nb gpu -pme gpu -bonded gpu -update gpu \
	-g ex3.4_${SLURM_NTASKS}x${OMP_NUM_THREADS}_${NPME}PME_jID${SLURM_JOB_ID} -nsteps -1 -maxh 0.017 -resethway -notunepme
   
```

:::warning
* Look at the absolute performance in the log file
* How does the absolute performance change when the number of PME ranks is changed? What is the optimal ration of PP to PME ranks for this system?
:::

## 4. Ensemble parallelization across multiple GPUs
==CB: This section seems very clear to me!==
:::success
:dart: **Learning goals**
* Understand how to set up and run ensemble simulations with `-multidir`.
* Understand the tradeoffs between simulation and aggregate ensemble throughput and how it relates to hardware utilization efficiency.
* _Advanced_: Explore the impact of task mapping on aggregate throughput.
:::


The `mdrun` simulation engine provides the _multi-simulation_ feature to run various types of ensemble simulations. For details of how to set up an ensemble simulation see [the GROMACS user guide section on multi-simulations](https://manual.gromacs.org/current/user-guide/mdrun-features.html#running-multi-simulations).

In this exercise we will learn how to run ensemble simulations using `mdrun` _multi-simulation_ feature. We will explore how to optimize ensemble performance and efficiency, and learn about tradeoffs between simulation throughput, aggregate throughput and how these relate to hardware utilization efficiency of heterogeneous hardware and GPU accelerators.

The example system uses a strongly coupled ensemble setup based on multi-walker AWH. The AWH setup is flexible in terms of ensemble size and can employ up to 32 members.

To run multi-simulations the mdrun option `-multidir` can be used. Note that the`-multidir` feature requires one input directory per ensemble member containing simulation input files (e.i tpr files) and this is where outputs will be written. You find the  directory structure in the [input tarball](https://a3s.fi/gmx-lumi/workshop-input-files.tar.gz) or in the folder /projappl/project_465000934/aquaporin/ensemble_inputs on LUMI.
==AV to check that the directory path is correct==


### Exercise 4.1: Ensemble runs with `-multidir`

When the simulation system is relatively small, it may not be able to fully saturate modern HPC GPUs. In such cases, we can achieve better hardware utilization by assigning multiple ensemble members to a GPU. By doing so we provide more (indenpendent) work to each GPU, which can significantly improve _aggregate simulation throughput_ and can allow making _more efficient use of GPU hardware_.

To explore this, we will use a fixed amount of hardware resources and vary the ensemble size. We will start with the GPU-resident setup from Exercise 2.2 and run multi-GPU ensemble runs on all eight GPUs of a single LUMI-G node.
```bash=
#!/bin/bash
#SBATCH --partition=small-g            # partition to use
#SBATCH --account=project_465000934    # project for billing
#SBATCH --reservation=gromacs_thursday # reservation for 25 January
#SBATCH --exclusive                    # new! to reserve the whole node
#SBATCH --time=00:10:00                # maximum execution time of 10 minutes
#SBATCH --nodes=1                      # we run on 1 node
#SBATCH --gpus-per-node=..             # fill in the number of GPU devices per node 
#SBATCH --ntasks-per-node=..           # fill in the number of  MPI rank

module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu
source ${GMXBIN}/lumi-affinity.sh   # 

export OMP_NUM_THREADS=             # fill in the value          
num_multi=....                      # change ensemble dimension 
#start with 8-way ensemble - update values above together with ?? in the command line below

srun --cpu-bind=${CPU_BIND} ./select_gpu gmx_mpi mdrun -multidir sim_{01..??} \
                   -nb gpu -pme gpu -bonded gpu -update gpu \
                   -g ex4.1_${SLURM_NNODES}N_multi${num_multi}_jID${SLURM_JOB_ID} \
                   -nsteps -1 -maxh 0.017 -resethway -notunepme
```
==AV to SP: please check that the above (CSC) is what you want. Below (your example). If yes delete the script below. I guess you want that the change also N of THREADS, good to give some hint==  
```bash=
#SBATCH --time=00:10:00
#SBATCH --nodes=1            
#SBATCH --ntasks=...         # Change this 
#SBATCH --cpus-per-task=...  # and this value as well, together with the "..." below

# module load ... 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

num_multi=16

mpirun -np $num_multi gmx_mpi mdrun \
   -multidir repl_{01..16} \
   -nb gpu -pme gpu -bonded gpu -update gpu \
   -g ex4.1_${SLURM_NNODES}N_multi${num_multi}_jID${SLURM_JOB_ID} \ 
   -nsteps -1 -maxh 0.017 -resethway -notunepme                    

```
* As a baseline, launch one simulation per GPU, hence an 8-way ensemble on one LUMI-G node. 
* Next, submit jobs with multiple (2,3 and 4) simulations per GPU on one LUMI-G node.

:::warning 
* How does the performance (ns/day) of _each ensemble member_ simulation change as you increase the number of simulations per GPU?
* How does the aggregate performance per node change as you increase the number of simulations per GPU?
* _Bonus_: log in to the compute node and observe the GPU utilization (using the `rocm-smi` tool) during the ensemble runs with lowest/highest aggregate performance.
:::



### Exercise 4.2: Trading efficiency for increasing simulation throughput 

In some cases, we want to increase the (non-aggregate) performance of the ensemble, e.g. if we need to sample longer. In such cases, we might have to trade performance for efficiency. 

In this exercise, we will combine what we have learnt in the exercise serie 3 and in ex. 4.1. We will scale each ensemble member across multiple GPUs to improve simulation performance. To do that, we will use a fixed 16-way ensemble and run it varying the amount of compute resources used.

```bash=
#!/bin/bash
#SBATCH --partition=small-g            # partition to use
#SBATCH --account=project_465000934    # project for billing
#SBATCH --reservation=gromacs_thursday # reservation for 25 January
#SBATCH --exclusive                    # new! to reserve the whole node
#SBATCH --time=00:10:00                # maximum execution time of 10 minutes
#SBATCH --nodes=...                    # fill in the number of node 
#SBATCH --gpus-per-node=..             # fill in the number of GPU devices per node 
#SBATCH --ntasks-per-node=..           # fill in the number of  MPI rank

module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu
source ${GMXBIN}/lumi-affinity.sh     # 

export OMP_NUM_THREADS=               # fill in the number of omp threads

export MPICH_GPU_SUPPORT_ENABLED=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1

num_multi=16                         # ensemble dimension 

srun --cpu-bind=${CPU_BIND} ./select_gpu \
     gmx_mpi mdrun -multidir sim_{01..16} \
                   -npme 1 \
                   -nb gpu -pme gpu -bonded gpu -update gpu \
                   -g ex4.2_${SLURM_NNODES}N_multi${num_multi}_jID${SLURM_JOB_ID} \
                   -nsteps -1 -maxh 0.017 -resethway -notunepme

```
==AV to SP in the original script in the command line was -npme ... (see below) while in the LUMI script (above) I set -npme 1 to avoid to many choice, please change back if you do not agree==
```bash=
#Command line in Szilard script
mpirun gmx_mpi mdrun \
   -multidir ${num_multi} repl_{01..16} \
   -npme ... \
   -nb gpu -pme gpu -bonded gpu -update gpu \
   -g ex4.2_${SLURM_NNODES}N_multi${num_multi}_jID${SLURM_JOB_ID} \ 
   -nsteps -1 -maxh 0.017 -resethway -notunepme                    
```

* As a baseline, we will use the single-node case, that is 16-way ensemble with 2 simulations per GPU (or a half-GPU per ensemble member) on a single LUMI-G node. Run this setup or you can reuse this result from the previous exercise.
* Next, submit jobs with increasing the amount of GPU resources assigned to each ensemble member, by increasing the total number of nodes used from 1 to 2, and 4.
* _Note_ that as mentioned earlier, this is a small simulation system and therefore its scaling is limited, hence it is not useful to try to use more than 2-4 GPUs/simulation.

:::warning 
* How does the performance (ns/day) of _each ensemble member_ simulation change as you increase the number of nodes/GPUs used?
* How does the aggregate performance per node change as you increase the number of nodes/GPUs used?
* _Bonus_: log in to the compute nodes and observe the GPU utilization using `rocm-smi` during the 1- and 4-node ensemble runs.
:::


### _Advanced Exercise 4.3_: Exploring task mapping for ensemble runs

In previous exercises, we used the fully GPU-offloaded GPU-resident mode, hence the CPU was left mostly idle except for the small amount of AWH and bonded computation which is not offloadable. As we explored in Exercise 2.2, exploiting the heterogeneous nature of the GROMACS engine, we can shift more work to the CPU cores that are otherwise left mostly idle.

Starting with the setup from the previous exercise, explore mapping some lighter-weight compute tasks to the CPU, e.g. `-bonded` or `-update`.

* As a baseline, we will use the fully GPU-offloaded runs from the previous exercise.
* Submit jobs with `-bonded cpu` and/or `-update cpu` 
* Bonus: repeat the same for a larger simulation system, like STMV . 
:::warning 
* Compare the log file when different tasks are GPU-offloaded. How do _per-simulation_ and _per-node_ aggregate performance change? 
* Bonus: Look at the STMV log files. Do a STMV and aquaporin show the same behavior?
::: 

[//]: <> (The following does not make sense on LUMI, so it's skipped)
[//]: <> (* Extra bonus: explore using multiple GPU runs mapped to the same set of GPUs e.g. 8-way ensemble per node, 2 GPUs/member, 2 members sharing a pair of GPUs. Try swapping around the order of PP and PME ranks to balance load on the shared GPUs.)
:::info
**License**
This material is shared under CC BY-SA 4.0 license. 
![CC BY-SA](https://hackmd.io/_uploads/SJvv38tuT.png "title" =240x84)


[//]: <> (This is a comment, it will not be included)
[//]: <> (in  the output file unless you use it in)
[//]: <> (a reference style link.)

[//]: # (This may be the most platform independent comment)

