#include "KorobovGenerator.hpp"
#include "VelocityGrid.hpp"
#include "CollisionNodes.hpp"
#include <cstdlib>



CollisionNodes::CollisionNodes(int n_nodes, VelocityGrid* p_velocity_grid) {
    this->n_nodes = n_nodes;
    this->p_velocity_grid = p_velocity_grid;
    collisions = std::vector<double>(8 * n_nodes);
    rel_velocities = std::vector<double>(3 * n_nodes);
    thetas = std::vector<double>(n_nodes);
    is_active = std::vector<bool>(n_nodes);
}
void CollisionNodes::randomizeNodes(std::vector<unsigned int> korobov_coefs) {
    for (int curr_node = 0; curr_node < n_nodes; curr_node++) {
        for (int i = 0; i < 8; i++) {
            // change to Korobov generated nodes
            collisions.push_back(std::rand() / static_cast<double>(RAND_MAX)); // uniform [0, 1]
        }
    }
}

void CollisionNodes::scaleValues(double v_cut, double s_max) {
    for (int curr_node = 0; curr_node < n_nodes; curr_node++) {
        for (int i = 0; i < 8; i++) {
            if (i < 6) {
                collisions[curr_node * 8 + i] = p_velocity_grid->gridLock(2 * v_cut * collisions[curr_node * 8 + i] - v_cut, i);
            }
            else if (i == 6) {
                collisions[curr_node * 8 + i] *= s_max;
                //got all necessary to calculate rel_velocities (g)
                rel_velocities[curr_node * 3] = collisions[curr_node * 8 + 3] - collisions[curr_node * 8];
                rel_velocities[curr_node * 3] = collisions[curr_node * 8 + 4] - collisions[curr_node * 8 + 1];
                rel_velocities[curr_node * 3] = collisions[curr_node * 8 + 5] - collisions[curr_node * 8 + 2];
                //got all necessary to calculate thetas
                thetas[curr_node] = acos(sqrt(collisions[curr_node * 8 + i]));
            }
            else {
                collisions[curr_node * 8 + i] *= 2 * M_PI;
            }
        }
    }
}

void CollisionNodes::checkOutOfSphere() {
    double v_c2 = pow(p_velocity_grid->getVCut(), 2);
    for (int curr_node = 0; curr_node < n_nodes; curr_node++) {
        double v_a2 = pow(collisions[curr_node], 2) + pow(collisions[curr_node + 1], 2) + pow(collisions[curr_node + 2], 2);
        double v_b2 = pow(collisions[curr_node + 3], 2) + pow(collisions[curr_node + 4], 2) + pow(collisions[curr_node + 5], 2);
        if ((v_a2 <= v_c2) && (v_b2 <= v_c2)) {
            is_active[curr_node] = true;        
        }
    }
}

void CollisionNodes::recalculateRelVelocities() {
    
}