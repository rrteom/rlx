#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include <fstream>

struct VectorVelocity {
    double x, y, z;
    VectorVelocity(double x, double y, double z) : x(x), y(y), z(z) {};
    VectorVelocity operator+(const VectorVelocity& b) const {
        return VectorVelocity(x + b.x, y + b.y, z + b.z);
    }; 
    VectorVelocity operator-(const VectorVelocity& b) const {
        return VectorVelocity(x - b.x, y - b.y, z - b.z);
    };
    VectorVelocity operator/(const double& d) const {
        return VectorVelocity(x / d, y / d, z / d);
    } 
    double pow2() const {
        return pow(x, 2) + pow(y, 2) + pow(z, 2);
    }
};

struct VectorIndex {
    int i, j, k;
    VectorIndex() {};
    VectorIndex(int i, int j, int k) : i(i), j(j), k(k) {};
    VectorIndex operator+(const VectorIndex& other) const {
        return VectorIndex(i + other.i, j + other.j, k + other.k);
    }
    VectorIndex operator-(const VectorIndex& other) const {
        return VectorIndex(i - other.i, j - other.j, k - other.k);
    }
    bool operator==(const VectorIndex& b) const {
        if ((i == b.i) && (j == b.j) && (k == b.k))
            return true;
        else
            return false;
    }; 
};

struct InterpNodes {
    double r;
    VectorIndex a, b, l, ls, m, ms;
    InterpNodes() {};
    InterpNodes(double r) : r(r) {};

};

double distrFunction(double x, double y, double z, double v_cut, double T, char type);

class VelocityGrid {
private:
    int Nx, Ny, Nz;
    double v_cut;
    std::vector<double> v_x, v_y, v_z; 
    std::vector<double> distr;

    inline int gridIndex(int i, int j, int k);
    inline int gridIndex(VectorIndex ix);
public:
    VelocityGrid(int N, double v_cut);
    VelocityGrid(int Nx, int Ny, int Nz, double v_cut);
    int countN0();
    void initDistr(double temp, char distr_type);
    void saveToFile(std::string path);
    double gridLockX(double v);
    double gridLockY(double v);
    double gridLockZ(double v);
    double gridLock(double v, int i);
    double getVCut();
    VectorVelocity gridLock(VectorVelocity vv);
    VectorIndex getClosestVeloctyIx(VectorVelocity vv);
    InterpNodes getInterpNodes(VectorVelocity v_a, VectorVelocity v_b, VectorVelocity v_a_new);
    VectorVelocity getVelocityByIx(VectorIndex ix);
    double calculateOmega(InterpNodes nd);
    inline double getDistr(VectorIndex ix);
};