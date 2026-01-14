#include <iostream>
#include <chrono>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include <unistd.h>   // gethostname
#include <cstdio>     // printf

#include "utils.h"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    auto start_parallel = std::chrono::high_resolution_clock::now();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // DEBUG MPI
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("Rank %d running on %s\n", rank, hostname);

    Grid grid = read_grid("../data/grid.pts");
    std::vector<Atom> ligand_atoms = read_xyz("../data/ligand.xyz");

    const int total_poses = 1000000;

    std::vector<std::vector<Atom>> poses;
    poses.reserve(total_poses);
    for (int i = 0; i < total_poses; ++i) {
        poses.push_back(transform_ligand(ligand_atoms, i));
    }

    double local_min = 1e9;
    int local_best = -1;

    #pragma omp parallel
    {
        double thread_min = 1e9;
        int thread_best = -1;
        
        // DEBUG OpenMP
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        #pragma omp critical
        printf("Rank %d: thread %d of %d\n", rank, tid + 1, nthreads);
        
        #pragma omp for schedule(static)
        for (int i = 0; i < poses.size(); i++) {
            const auto& pose = poses[i];

            double total = 0.0;
            for (const auto& atom : pose) {
                total += trilinear_interp(grid, atom.x, atom.y, atom.z);
            }

            if (total < thread_min) {
                thread_min = total;
                thread_best = i;
            }
        }

        #pragma omp critical
        {
            if (thread_min < local_min) {
                local_min = thread_min;
                local_best = thread_best;
            }
        }
    }

    auto end_parallel = std::chrono::high_resolution_clock::now();

    if (rank == 0) {
        std::chrono::duration<double> elapsed = end_parallel - start_parallel;
        std::cout << "Interpolated value = " << local_min << std::endl;
        std::cout << "Best pose = " << local_best << std::endl;
        std::cout << "Time taken = " << elapsed.count() << std::endl;
    }
    MPI_Finalize();
    return 0;
}