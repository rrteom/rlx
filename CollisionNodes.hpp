#include <vector>
#include <utility>
#include <tuple>
#include "VelocityGrid.hpp"

class CollisionNodes {
private:
    int n_nodes;
    std::vector<double> collisions, rel_velocities, thetas;
    std::vector<bool> is_active;
    VelocityGrid* p_velocity_grid;
    bool g_recalculated;
public:
    CollisionNodes(int n_nodes, VelocityGrid* p_velocity_grid);
    void randomizeNodes(std::vector<unsigned int> korobov_coefs);
    void scaleValues(double v_cut, double S_max);
    void checkOutOfSphere();
    void recalculateRelVelocities();
    std::pair<VectorVelocity, VectorVelocity> afterCollisionVelocities(int collision_no); // just for testing
    std::tuple<VectorVelocity, VectorVelocity, VectorVelocity> getVelocitiesForInterp(int coll_no); 
};