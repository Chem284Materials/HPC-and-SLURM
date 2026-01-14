# 💻 CHEM 284 - HPC and SLURM

## 🧪 Goal

The goal of this lab is to:

1. Learn about **high performace computing**.
2. Learn how to use **SLURM to run jobs on a HPC cluster**. 
3. Practice using **a hybrid parallel scheme on a molecular docking virtual screen**.

---
## 🗂️ Provided

- A `docker-compose` file to set up the HPC cluster environment.
- A `CMakeLists.txt` file to help build the executables.
- A `data` directory with the relevant files.
- A `src` directory with the relevant files.

---
## 💻 Setup
```bash
./build_images.sh # You may need to chmod +x
docker-compose up -d

Creating network "lab_5_default" with the default driver
Creating lab_5_slurmjupyter_1 ... done
Creating lab_5_slurmmaster_1  ... done
Creating lab_5_slurmnode2_1   ... done
Creating lab_5_slurmnode3_1   ... done
Creating lab_5_slurmnode1_1   ... done

# Navigate to
http://localhost:8888/lab?

# To finish
docker-compose down
```

This lab is a little different and instead of working directly inside a container you will be creating a bunch of docker containers stitched together to form a cluster using `docker-compose`. The main entry point is the `slurm-jupyter-hpc` container which we access using `http://localhost:8888/lab?`. Once you navigate to the jupyter lab page, you will find a shared folder with all the necessary files for the lab. You can create a terminal session by clicking File -> New -> Terminal. This terminal is inside the `slurm-jupyter-hpc` container and will allow you run all your commands and build using CMakeLists.txt. Additionally, you can navigate to any file in the `shared` directory and open it and edit it.

### Validate setup

Once you have a terminal session go into the `shared` directory and run the job.sh script under SLURM, using `sbatch job.sh`. It should submit and a few seconds later you should see a result.out file that should show you all 3 of your nodes working with 2 tasks per node. Look at `job.sh` for more info on the format of the SLURM commands.

```bash
admin@slurmjupyter:~/shared$ sbatch job.sh 
Submitted batch job 8

admin@slurmjupyter:~/shared$ cat result.out 
Process started 2026-01-13 05:02:26.875649
NODE : slurmnode1
PID  : 769
CORES: 16
Executing for 15 secs
Process finished 2026-01-13 05:02:41.890953

Process started 2026-01-13 05:02:26.875778
NODE : slurmnode1
PID  : 770
CORES: 16
Executing for 15 secs
Process finished 2026-01-13 05:02:41.891062

Process started 2026-01-13 05:02:26.899233
NODE : slurmnode2
PID  : 474
CORES: 16
Executing for 15 secs
Process finished 2026-01-13 05:02:41.913445

Process started 2026-01-13 05:02:26.899232
NODE : slurmnode3
PID  : 277
CORES: 16
Executing for 15 secs
Process finished 2026-01-13 05:02:41.913396

Process started 2026-01-13 05:02:26.900015
NODE : slurmnode2
PID  : 473
CORES: 16
Executing for 15 secs
Process finished 2026-01-13 05:02:41.915644

Process started 2026-01-13 05:02:26.900064
NODE : slurmnode3
PID  : 278
CORES: 16
Executing for 15 secs
Process finished 2026-01-13 05:02:41.915731
```

## ✅ Tasks
### Hybrid MPI-OpenMP docking
So far we have parallelized a virtual screen with MPI using a receptor, and with OpenMP using a grid. Now that we have a cluster we will try a hybrid approach using MPI across nodes and OpenMP within them. The `shared/src/hybrid_docking.cpp` file already has the OpenMP implementation done for you and your goal is to add the MPI components.

1) To keep it simple each node will read in the grid and thus we will be duplicate the grid across our nodes.
2) Since we are not reading in poses but creating them with `transform_ligand()` each node can start with the same starting ligand and then when we divide the work up for each rank, we can get start and end indices for the number of poses each rank is responsible for. Then on each node we can just generate their local poses with the indices. Make sure to account for the case when you have an even number of poses and an odd number of nodes.
3) You will need a synchronization step so that all nodes have their respective poses created before you move on to the OpenMP stage.
4) Once the OpenMP section is done you will need to convert the local best pose index for each node to a global best pose index across all nodes.
5) The last step is to do a final reduction on the global best pose index across all nodes.

Build as you normally do in the jupyter-lab terminal using cmake:
```bash
mkdir build
cd build
cmake ..
make
```

To run the code you have 2 options and you should try both of them:

1) Using `salloc` to get an allocation from SLURM to run your code. You can then simply call `mpirun ./hybrid_docking` and it will run it on your allocated nodes. `salloc` uses the same SLURM keywords as `sbatch`: 

`--nodes`: Number of nodes to run on
`--ntasks-per-node`: Number of MPI ranks to run on each node
`--cpus-per-task`: CPU reservation, Number of OpenMP threads

Note: You will need to set OMP_NUM_THREADS before you run with mpirun since by default even if we asked for more than 1 cpu-per-task, it will default to 1 on our cluster or be empty.
```bash
admin@slurmjupyter:~/shared/build$ salloc --nodes=2 --ntasks-per-node=1 --cpus-per-task=4
salloc: Granted job allocation 9

admin@slurmjupyter:~/shared/build$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
                 9  slurmpar     bash    admin  R       0:02      2 slurmnode[1-2]

admin@slurmjupyter:~/shared/build$ echo $OMP_NUM_THREADS
1

admin@slurmjupyter:~/shared/build$ mpirun --bind-to none --map-by slot ./hybrid_docking
Rank 1 running on slurmnode2
Rank 0 running on slurmnode1
Rank: 1 count: 500000
Rank: 0 count: 500000
Rank 1: thread 1 of 1
Rank 0: thread 1 of 1
Interpolated value = -6.38924
Best pose = 306172
Time taken = 0.719306
```

2) Using `sbatch` you can create a bash script using the `#SBATCH` directives to define the nodes, ntasks-per-node and cpus-per-task as well as set an output file for the stdout to log to. In this script you will also use mpirun to run your executable.

Note: When using sbatch with our code you will need to have the script inside of the build directory to run from.

### Extra time
Try different nodes, cpus-per-task configurations, and a larger or smaller number of poses. How does it scale? Are there any limitations you see?

#### References
Link to medium article explaining how to set up a docker-compose based [HPC cluster](https://medium.com/analytics-vidhya/slurm-cluster-with-docker-9f242deee601).
