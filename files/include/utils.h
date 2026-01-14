#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

struct Atom {
    double x, y, z;
};

std::vector<Atom> read_xyz(std::string file_name) {
    std::ifstream infile(file_name);
    
    std::vector<Atom> atoms;
    if (!infile) {
        std::cerr << "Error: could not open file, atoms will be empty\n";
        return atoms;
    }

    // First line is number of atoms
    int num_atoms;
    infile >> num_atoms;

    double x, y, z;
    while(infile >> x >> y >> z) {
        atoms.push_back({x, y, z});
    }
    return atoms;
}

inline void euler_to_matrix(double alpha, double beta, double gamma,
                            double R[3][3])
{
    double ca = cos(alpha), sa = sin(alpha);
    double cb = cos(beta),  sb = sin(beta);
    double cg = cos(gamma), sg = sin(gamma);

    R[0][0] =  cb*cg;
    R[0][1] = -cb*sg;
    R[0][2] =  sb;

    R[1][0] =  sa*sb*cg + ca*sg;
    R[1][1] = -sa*sb*sg + ca*cg;
    R[1][2] = -sa*cb;

    R[2][0] = -ca*sb*cg + sa*sg;
    R[2][1] =  ca*sb*sg + sa*cg;
    R[2][2] =  ca*cb;
}

std::vector<Atom> transform_ligand(const std::vector<Atom>& ligand,
                                   int pose_index)
{
    std::vector<Atom> out = ligand;

    // -----------------------------------------------------------
    // 1) Compute ligand centroid
    // -----------------------------------------------------------
    double cx = 0, cy = 0, cz = 0;
    for (auto& a : ligand) {
        cx += a.x;
        cy += a.y;
        cz += a.z;
    }
    cx /= ligand.size();
    cy /= ligand.size();
    cz /= ligand.size();

    // -----------------------------------------------------------
    // 2) Generate Euler rotation matrix
    // -----------------------------------------------------------
    double alpha = (pose_index * 0.37);
    double beta  = (pose_index * 0.51);
    double gamma = (pose_index * 0.29);

    double R[3][3];
    euler_to_matrix(alpha, beta, gamma, R);

    // -----------------------------------------------------------
    // 3) Translation
    // -----------------------------------------------------------
    double tx = 5.0 * sin(pose_index * 0.21);
    double ty = 5.0 * cos(pose_index * 0.13);
    double tz = 5.0 * sin(pose_index * 0.17 + 0.5);

    // -----------------------------------------------------------
    // 4) Apply:   (p - centroid) → rotate → + centroid → + translation
    // -----------------------------------------------------------
    for (auto& a : out) {
        double x = a.x - cx;
        double y = a.y - cy;
        double z = a.z - cz;

        double rx = R[0][0]*x + R[0][1]*y + R[0][2]*z;
        double ry = R[1][0]*x + R[1][1]*y + R[1][2]*z;
        double rz = R[2][0]*x + R[2][1]*y + R[2][2]*z;

        // Re-center and translate
        a.x = rx + cx + tx;
        a.y = ry + cy + ty;
        a.z = rz + cz + tz;
    }

    return out;
}

/*
 File format:

 n
 ligand_center_x ligand_center_y ligand_center_z
 x_min y_min z_min
 dx
 n^3
 ix iy iz value
 ...
*/

struct Grid {
    int n;                     // points per dimension
    double x_min, y_min, z_min;
    double dx;
    std::vector<double> values; // flattened grid, size = n^3
};  

Grid read_grid(const std::string& filename)
{
    std::ifstream f(filename);
    if (!f) {
        throw std::runtime_error("Could not open grid file");
    }

    Grid g;
    double cx, cy, cz;
    int total_points;

    f >> g.n;
    f >> cx >> cy >> cz;                    // ligand center (not needed later)
    f >> g.x_min >> g.y_min >> g.z_min;
    f >> g.dx;
    f >> total_points;

    if (total_points != g.n * g.n * g.n) {
        throw std::runtime_error("Grid size mismatch");
    }

    g.values.resize(total_points);

    int ix, iy, iz;
    double val;
    while (f >> ix >> iy >> iz >> val) {
        int idx = ix + g.n * (iy + g.n * iz);
        g.values[idx] = val;
    }

    return g;
}

/* ------------------- Flattened indexing ------------------- */
inline double grid_value(const Grid& g, int ix, int iy, int iz)
{
    return g.values[ix + g.n * (iy + g.n * iz)];
}

double trilinear_interp(const Grid& g,
                        double x, double y, double z)
{
    // Convert coordinates to floating-point indices
    double i_f = (x - g.x_min) / g.dx;
    double j_f = (y - g.y_min) / g.dx;
    double k_f = (z - g.z_min) / g.dx;

    int i0 = std::floor(i_f);
    int j0 = std::floor(j_f);
    int k0 = std::floor(k_f);

    // Clamp so we never go out of bounds
    i0 = std::max(0, std::min(i0, g.n - 2));
    j0 = std::max(0, std::min(j0, g.n - 2));
    k0 = std::max(0, std::min(k0, g.n - 2));

    int i1 = i0 + 1;
    int j1 = j0 + 1;
    int k1 = k0 + 1;

    double xd = i_f - i0;
    double yd = j_f - j0;
    double zd = k_f - k0;

    // Corner values
    double C000 = grid_value(g, i0, j0, k0);
    double C100 = grid_value(g, i1, j0, k0);
    double C010 = grid_value(g, i0, j1, k0);
    double C110 = grid_value(g, i1, j1, k0);
    double C001 = grid_value(g, i0, j0, k1);
    double C101 = grid_value(g, i1, j0, k1);
    double C011 = grid_value(g, i0, j1, k1);
    double C111 = grid_value(g, i1, j1, k1);

    // Interpolate along x
    double C00 = C000 * (1 - xd) + C100 * xd;
    double C01 = C001 * (1 - xd) + C101 * xd;
    double C10 = C010 * (1 - xd) + C110 * xd;
    double C11 = C011 * (1 - xd) + C111 * xd;

    // Interpolate along y
    double C0 = C00 * (1 - yd) + C10 * yd;
    double C1 = C01 * (1 - yd) + C11 * yd;

    // Interpolate along z
    return C0 * (1 - zd) + C1 * zd;
}