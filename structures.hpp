#ifndef STRUCTURES_HEADER
#define STRUCTURES_HEADER

#include <cmath>

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

#endif