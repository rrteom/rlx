#include "KorobovGenerator.hpp"
#include "VelocityGrid.hpp"
#include "CollisionNodes.hpp"
#include <cstdlib>

#include <iostream>



CollisionNodes::CollisionNodes(int n_nodes, VelocityGrid* p_velocity_grid, double s_max) {
    this->n_nodes = n_nodes;
    this->p_velocity_grid = p_velocity_grid;
    this->s_max = s_max;
    collisions = std::vector<double>(8 * n_nodes);
    rel_velocities = std::vector<double>(3 * n_nodes);
    thetas = std::vector<double>(n_nodes);
    is_active = std::vector<bool>(n_nodes);
    g_recalculated = false;
}
void CollisionNodes::randomizeNodes(std::vector<unsigned int> korobov_coefs) {
    for (int curr_node = 0; curr_node < n_nodes; curr_node++) {
        for (int i = 0; i < 8; i++) {
            // change to Korobov generated nodes
            collisions[curr_node * 8 + i] = std::rand() / static_cast<double>(RAND_MAX); // uniform [0, 1]
        }
    }
}

void CollisionNodes::scaleValues() {
    double v_cut = p_velocity_grid->getVCut();
    for (int curr_node = 0; curr_node < n_nodes; curr_node++) {
        for (int i = 0; i < 8; i++) {
            if (i < 6) {
                collisions[curr_node * 8 + i] = p_velocity_grid->gridLock(2 * v_cut * collisions[curr_node * 8 + i] - v_cut, i);
            }
            else if (i == 6) {
                collisions[curr_node * 8 + i] *= s_max;
                //got all necessary to calculate rel_velocities (g)
                rel_velocities[curr_node * 3] = collisions[curr_node * 8 + 3] - collisions[curr_node * 8];
                rel_velocities[curr_node * 3 + 1] = collisions[curr_node * 8 + 4] - collisions[curr_node * 8 + 1];
                rel_velocities[curr_node * 3 + 2] = collisions[curr_node * 8 + 5] - collisions[curr_node * 8 + 2];
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
        double v_a2 = pow(collisions[8 * curr_node], 2) + pow(collisions[8 * curr_node + 1], 2) + pow(collisions[8 * curr_node + 2], 2);
        double v_b2 = pow(collisions[8 * curr_node + 3], 2) + pow(collisions[8 * curr_node + 4], 2) + pow(collisions[8 * curr_node + 5], 2);
        if ((v_a2 <= v_c2) && (v_b2 <= v_c2)) {
            is_active[curr_node] = true;      
        }
    }
}

void CollisionNodes::recalculateRelVelocities() {
    for (int curr_node = 0; curr_node < n_nodes; curr_node++) {
        double g = sqrt(pow(rel_velocities[3 * curr_node], 2) + pow(rel_velocities[3 * curr_node + 1], 2) + pow(rel_velocities[3 * curr_node + 2], 2));
        double g_xy = sqrt(pow(rel_velocities[3 * curr_node], 2) + pow(rel_velocities[3 * curr_node + 1], 2));
        double sin_e = sin(collisions[curr_node * 8 + 7]), cos_e = cos(collisions[curr_node * 8 + 7]); 
        double sin_th = sin(thetas[curr_node]), cos_th = cos(thetas[curr_node]);
        double g_x = rel_velocities[3 * curr_node], g_y = rel_velocities[3 * curr_node + 1], g_z = rel_velocities[3 * curr_node + 2];
        if (g_xy == 0) {
            double new_x = g_x * sin_e * sin_th;
            double new_y = g_y * cos_e * sin_th;
            double new_z = g_z * cos_th;
            rel_velocities[3 * curr_node] = new_x;
            rel_velocities[3 * curr_node + 1] = new_y;
            rel_velocities[3 * curr_node + 2] = new_z; 
        }
        else {
            double new_x = g_x * cos_th + sin_th / g_xy * (-g_x * g_z * cos_e + g * g_y * sin_e);
            double new_y = g_y * cos_th + sin_th / g_xy * (-g_y * g_z * cos_e - g * g_x * sin_e);
            double new_z = g_z * cos_th + g_xy * cos_e * sin_th;
            rel_velocities[3 * curr_node] = new_x;
            rel_velocities[3 * curr_node + 1] = new_y;
            rel_velocities[3 * curr_node + 2] = new_z; 
        }
    }
    g_recalculated = true;
    return;
}

std::pair<VectorVelocity, VectorVelocity> CollisionNodes::afterCollisionVelocities(int coll_no) {
    // just for testing
    if (!g_recalculated) {
        recalculateRelVelocities();
    }
    VectorVelocity first(
        (collisions[coll_no * 8] + collisions[coll_no * 8 + 3] - rel_velocities[coll_no * 3]) / 2,
        (collisions[coll_no * 8 + 1] + collisions[coll_no * 8 + 4] - rel_velocities[coll_no * 3 + 1]) / 2,
        (collisions[coll_no * 8 + 2] + collisions[coll_no * 8 + 5] - rel_velocities[coll_no * 3 + 2]) / 2
    );
    VectorVelocity second(
        (collisions[coll_no * 8] + collisions[coll_no * 8 + 3] + rel_velocities[coll_no * 3]) / 2,
        (collisions[coll_no * 8 + 1] + collisions[coll_no * 8 + 4] + rel_velocities[coll_no * 3 + 1]) / 2,
        (collisions[coll_no * 8 + 2] + collisions[coll_no * 8 + 5] + rel_velocities[coll_no * 3 + 2]) / 2
    );
    return std::make_pair(first, second);
}

std::tuple<VectorVelocity, VectorVelocity, VectorVelocity> CollisionNodes::getVelocitiesForInterp(int coll_no) {
    VectorVelocity v_a(collisions[coll_no * 8], collisions[coll_no * 8 + 1], collisions[coll_no * 8 + 2]);
    VectorVelocity v_b(collisions[coll_no * 8 + 3], collisions[coll_no * 8 + 4], collisions[coll_no * 8 + 5]);
    VectorVelocity v_a_new(
        (collisions[coll_no * 8] + collisions[coll_no * 8 + 3] - rel_velocities[coll_no * 3]) / 2,
        (collisions[coll_no * 8 + 1] + collisions[coll_no * 8 + 4] - rel_velocities[coll_no * 3 + 1]) / 2,
        (collisions[coll_no * 8 + 2] + collisions[coll_no * 8 + 5] - rel_velocities[coll_no * 3 + 2]) / 2
    );
    return std::tuple(v_a, v_b, v_a_new);
}

bool CollisionNodes::isActive(int coll_no) {
    return is_active[coll_no];
}

void CollisionNodes::saveToFile(std::string path) {
    std::ofstream out_stream;
    out_stream.open(path);
    out_stream << n_nodes << ' ' << g_recalculated << std::endl;
    
    for (auto i: collisions) 
        out_stream << i << ' ';
    out_stream << std::endl;

    for (auto i: rel_velocities) 
        out_stream << i << ' ';
    out_stream << std::endl;

    for (auto i: thetas) 
        out_stream << i << ' ';
    out_stream << std::endl;

    for (auto i: is_active) 
        out_stream << i << ' ';
    out_stream << std::endl;

    out_stream.close();
    return;
}