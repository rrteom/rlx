#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include <fstream>
#include "VelocityGrid.hpp"


double distrFunction(double x, double y, double z, double v_cut, double temp, char distr_type) {
    double v_squared = pow(x, 2) + pow(y, 2) + pow(z, 2);
    if (distr_type == '1') 
        return v_squared <= pow(v_cut, 2) ? exp(-0.5 * v_squared) + pow(temp, -1.5) * exp(-0.5 * v_squared / temp) : 0.0; 
    else if (distr_type == '2') {
        if (v_squared > pow(v_cut, 2))
            return 0.0;
        else if (x > 0)
            return exp(-0.5 * v_squared);
        else
            return pow(temp, -1.5) * exp(-0.5 * v_squared / temp);
    }
    else 
        return v_squared <= pow(v_cut, 2) ? 1.0 : 0.0; // for test
}

VelocityGrid::VelocityGrid(int Nx, int Ny, int Nz, double v_cut) {
    this -> Nx = Nx;
    this -> Ny = Ny;
    this -> Nz = Nz;
    this -> v_cut = v_cut;
    for (int i = 0; i < Nx; i++) v_x.push_back(-v_cut + (i + 0.5) * 2 * v_cut / Nx);
    for (int i = 0; i < Ny; i++) v_y.push_back(-v_cut + (i + 0.5) * 2 * v_cut / Ny);
    for (int i = 0; i < Nz; i++) v_z.push_back(-v_cut + (i + 0.5) * 2 * v_cut / Nz);
    distr.resize(Nx * Ny * Nz);
}

VelocityGrid::VelocityGrid(int N, double v_cut) : VelocityGrid::VelocityGrid(N, N, N, v_cut) {}

int VelocityGrid::countN0() {
    // плохо, должно быть красивое математическое решение
    int count = 0;
    for (int i = 0; i < v_x.size(); i++) {
        for (int j = 0; j < v_y.size(); j++) {
            for (int k = 0; k < v_z.size(); k++) {
                if (pow(v_x[i], 2) + pow(v_y[j], 2) + pow(v_z[k], 2) <= pow(v_cut, 2))
                    count++;
            }
        }
    }
    return count;
}

inline int VelocityGrid::gridIndex(int x, int y, int z) {
    return x * Ny * Nz + y * Nz + z; 
} 

inline int VelocityGrid::gridIndex(VectorIndex ix) {
    return gridIndex(ix.i, ix.j, ix.k);
}

inline double VelocityGrid::getDistr(VectorIndex ix) {
    return distr[gridIndex(ix)];
}

void VelocityGrid::initDistr(double temp = 0.5, char distr_type = '1') {
    // init unnormalised
    double sum = 0;
    for (int i = 0; i < v_x.size(); i++) {
        for (int j = 0; j < v_y.size(); j++) {
            for (int k = 0; k < v_z.size(); k++) {
                double current_distr_value = distrFunction(v_x[i], v_y[j], v_z[k], v_cut, temp, distr_type);
                distr[gridIndex(i, j, k)] = current_distr_value;
                sum += current_distr_value;
            }
        }
    }
    // normalisation
    sum *= pow(2 * v_cut, 3) / (Nx - 1) / (Ny - 1) / (Nz - 1); 
    for (int i = 0; i < v_x.size(); i++) {
        for (int j = 0; j < v_y.size(); j++) {
            for (int k = 0; k < v_z.size(); k++) {
                if (distr[gridIndex(i, j, k)] != 0)
                    distr[gridIndex(i, j, k)] = distr[gridIndex(i, j, k)] / sum;
            }
        }
    }
}

void VelocityGrid::saveToFile(std::string path) {
    std::ofstream out_stream;
    out_stream.open(path);
    out_stream << Nx << ' ' << Ny << ' ' << Nz << ' ' << v_cut << std::endl;
    for (int i = 0; i < distr.size(); i++)
        out_stream << distr[i] << ' ';
    out_stream.close();
    return;
}

double VelocityGrid::gridLockY(double v) {
    int index = static_cast<int>(floor(Ny / 2 * (1 + v / v_cut)));
    if (index < 0) 
        return v_y[0];
    else if (index >= Ny) 
        return v_y[Ny - 1];
    else 
        return v_y[index];
}

double VelocityGrid::gridLockX(double v) {
    int index = static_cast<int>(floor(Nx / 2 * (1 + v / v_cut)));
    if (index < 0) 
        return v_x[0];
    else if (index >= Nx) 
        return v_x[Nx - 1];
    else 
        return v_x[index];
}

double VelocityGrid::gridLockZ(double v) {
    int index = static_cast<int>(floor(Nz / 2 * (1 + v / v_cut)));
    if (index < 0) 
        return v_z[0];
    else if (index >= Nz) 
        return v_z[Nz - 1];
    else 
        return v_z[index];
}

double VelocityGrid::gridLock(double v, int i) {
    if ((i == 0) || (i == 3))
        return gridLockX(v);
    else if ((i == 1) || (i == 4))
        return gridLockY(v);
    else if ((i == 2) || (i == 5))
        return gridLockZ(v);
    else 
        return v;
}

VectorVelocity VelocityGrid::gridLock(VectorVelocity vv) {
    return VectorVelocity(gridLockX(vv.x), gridLockY(vv.y), gridLockZ(vv.z));
}

VectorIndex VelocityGrid::getClosestVeloctyIx(VectorVelocity vv) {
    int i = static_cast<int>(floor(Nx / 2 * (1 + vv.x / v_cut)));
    int j = static_cast<int>(floor(Ny / 2 * (1 + vv.y / v_cut)));
    int k = static_cast<int>(floor(Nz / 2 * (1 + vv.z / v_cut)));
    return VectorIndex(i, j, k);
}

double VelocityGrid::getVCut() {
    return v_cut;
}

VectorVelocity VelocityGrid::getVelocityByIx(VectorIndex ix) {
    double x, y, z;
    if (ix.i < 0)
        x = -v_cut;
    else if (ix.i >= Nx)
        x = v_cut;
    else 
        x = v_x[ix.i];
    
    if (ix.j < 0)
        y = -v_cut;
    else if (ix.j >= Ny)
        y = v_cut;
    else 
        y = v_y[ix.j];

    if (ix.k < 0)
        z = -v_cut;
    else if (ix.k >= Nz)
        z = v_cut;
    else 
        z = v_z[ix.k];
    return VectorVelocity(x, y, z);
}

InterpNodes VelocityGrid::getInterpNodes(VectorVelocity v_a, VectorVelocity v_b, VectorVelocity v_a_new) {
    VectorVelocity v_near = gridLock(v_a_new), v_near_1 = gridLock(v_a + v_b - v_a_new);
    if ((v_near.pow2() > pow(v_cut, 2)) || (v_near_1.pow2() > pow(v_cut, 2))) {
        return InterpNodes(-1); // if r < 0 we deactivate (do not count) collision
    }
    VectorIndex eta_near = getClosestVeloctyIx(v_near), eta_near_1 = getClosestVeloctyIx(v_near_1);
    
    int i_center = static_cast<int>(floor(Nx / 2 * (1 + v_a_new.x / v_cut) - 0.5));
    int j_center = static_cast<int>(floor(Ny / 2 * (1 + v_a_new.y / v_cut) - 0.5));
    int k_center = static_cast<int>(floor(Nz / 2 * (1 + v_a_new.z / v_cut) - 0.5));
    VectorIndex eta_center(i_center, j_center, k_center);
    
    VectorVelocity v_cm = (v_a + v_b) / 2; // mass center velocity
    int q_near;
    std::vector<VectorIndex> eta_around;
    std::vector<double> energ_around;
    for (int x = 0; x <=1; x++) {
        for (int y = 0; y <=1; y++) {
            for (int z = 0; z <=1; z++) {
                VectorIndex eta_q(eta_center.i + x, eta_center.j + y, eta_center.k + z);
                eta_around.push_back(eta_q);
                energ_around.push_back((getVelocityByIx(eta_q) - v_cm).pow2());
                if (eta_q == eta_near) {
                    q_near = 4 * x + 2 * y + z;
                }
            }
        }
    }
    
    double energ_0 = (v_a - v_b).pow2() / 4;
    double energ_1, energ_2;
    InterpNodes result;
    result.a = getClosestVeloctyIx(v_a);
    result.b = getClosestVeloctyIx(v_b);

    if (fabs(energ_around[q_near] - energ_0) < 1e-18) {
        result.r = 1;
        result.l = eta_near;
        result.ls = eta_near;
        result.m = eta_near_1;
        result.ms = eta_near_1;
        return result;
    }
    else if (energ_around[q_near] > energ_0) {
        double min_distance_2 = -1;
        for (int q; q < eta_around.size(); q++) {
            if ((energ_around[q] > energ_0) || (getVelocityByIx(eta_around[q]).pow2() > pow(v_cut, 2)))
                continue;
            VectorIndex this_index = eta_around[q];
            VectorIndex pair_index = result.a + result.b - this_index;
            if (getVelocityByIx(pair_index).pow2() > pow(v_cut, 2))
                continue;
            double distance_2 = (v_a_new - getVelocityByIx(this_index)).pow2();
            if ((min_distance_2 < 0) || (distance_2 < min_distance_2)) {
                min_distance_2 = distance_2;
                result.l = this_index;
                result.m = pair_index;
                energ_1 = energ_around[q];
            }
        }
        if (min_distance_2 < 0) 
            return InterpNodes(-1); 
        result.ls = eta_near;
        result.ms = eta_near_1;
        energ_2 = energ_around[q_near];
    }
    else {
        double min_distance_2 = -1;
        for (int q; q < eta_around.size(); q++) {
            if ((energ_around[q] < energ_0) || (getVelocityByIx(eta_around[q]).pow2() > pow(v_cut, 2)))
                continue;
            VectorIndex this_index = eta_around[q];
            VectorIndex pair_index = result.a + result.b - this_index;
            if (getVelocityByIx(pair_index).pow2() > pow(v_cut, 2))
                continue;
            double distance_2 = (v_a_new - getVelocityByIx(this_index)).pow2();
            if ((min_distance_2 < 0) || (distance_2 < min_distance_2)) {
                min_distance_2 = distance_2;
                result.ls = this_index;
                result.ms = pair_index;
                energ_2 = energ_around[q];
            }
        }
        if (min_distance_2 < 0) 
            return InterpNodes(-1);
        result.l = eta_near;
        result.m = eta_near_1;
        energ_1 = energ_around[q_near];
    }
    result.r = (energ_0 - energ_1) / (energ_2 - energ_1);
    return result;
}

double VelocityGrid::calculateOmega(InterpNodes nodes) {
    double f_l_f_m = getDistr(nodes.l) * getDistr(nodes.m);
    double v_diff = sqrt((getVelocityByIx(nodes.a) - getVelocityByIx(nodes.b)).pow2());
    if (f_l_f_m == 0)
        return 0;
    return f_l_f_m * pow(getDistr(nodes.ls) * getDistr(nodes.ms) / f_l_f_m, nodes.r) * v_diff;
}