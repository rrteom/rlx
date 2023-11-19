#ifndef VELOCITY_GRID
#define VELOCITY_GRID

#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include <fstream>
#include "structures.hpp"

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
    InterpNodes getInterpNodes(std::tuple<VectorVelocity, VectorVelocity, VectorVelocity>);
    VectorVelocity getVelocityByIx(VectorIndex ix);
    double calculateOmega(InterpNodes nd);
    inline double getDistr(VectorIndex ix);
    void updateDistr(InterpNodes nodes, double omega, double constant);
    
    double getConcentration();
    double getEnergy();
    VectorVelocity getMomentum();
    void printIntegrals();

};

#endif