#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include <fstream>

static double distrFunction(double x, double y, double z, double T);

class VelocityGrid {
private:
    int Nx, Ny, Nz;
    double v_cut;
    std::vector<double> v_x, v_y, v_z; 
    std::vector<double> distr;

    inline int gridIndex(int i, int j, int k);
public:
    VelocityGrid(int N, double v_cut);
    VelocityGrid(int Nx, int Ny, int Nz, double v_cut);
    int countN0();
    void initDistr(double temp, char distr_type);
    void saveToFile(std::string path);
};

static double distrFunction(double x, double y, double z, double v_cut, double temp = 0.5, char distr_type = '1') {
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