#include <iostream>
#include <chrono>
#include <vector>

#include "utils.h"

int main()
{
    Grid grid = read_grid("../data/grid.pts");

    std::vector<Atom> ligand_atoms = read_xyz("../data/ligand.xyz");
    // Generate poses
    std::vector<std::vector<Atom>> poses;
    for (int i = 0; i < 1000000; i++) {
        poses.push_back(transform_ligand(ligand_atoms, i));
    }

    // Start timing
    auto start_serial = std::chrono::high_resolution_clock::now();
    // Find global minimum energy and pose
    double global_min = 10000.0;
    int best_pose = -1;
    for (int i = 0; i < poses.size(); i++) {
        double total = 0.0;
        for (auto& atom: poses[i]) {
            total += trilinear_interp(grid, atom.x, atom.y, atom.z);
        }
        if (total < global_min) {
            global_min = total;
            best_pose = i;
        }
    }

    // End timing
    auto end_serial = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_serial = end_serial - start_serial;
    double t_serial = elapsed_serial.count();

    std::cout << "Interpolated value = " << global_min << std::endl;
    std::cout << "Best pose = " << best_pose << std::endl;
    std::cout << "Time taken = " << t_serial << std::endl;

    return 0;
}