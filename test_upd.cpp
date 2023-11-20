#include "CollisionNodes.hpp"
#include "UtilsRandom.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

int main() {
    // params block
    double tau = 0.02, t_end = 2; // time interval and step
    double s_max = 1, v_cut=4.8, temperature = 0.95;
    int n_velocity_nodes = 20; // = Nx = Ny = Nz
    int n_coll_nodes = 200003; //N_{\nu}
    char distr_type = '2';
    srand(42); // needed for reproducible permutation generation
    std::vector<unsigned int> korobov_coefs = {1, 47369, 188507, 54145, 156036, 158419, 37051, 42494};
    // params end
    int active_collisions = 0, valid_collisions = 0;

    n_coll_nodes = next_prime(n_coll_nodes);

    VelocityGrid velocity_grid(n_velocity_nodes, v_cut);
    velocity_grid.initDistr(temperature, distr_type);
    velocity_grid.printIntegrals();
    velocity_grid.saveToFile("checkpoints_new/log_init" + std::to_string(distr_type) + ".txt");
    
    int n_0 = velocity_grid.countN0();
    double integr_constant = s_max * pow(2 * v_cut / n_velocity_nodes, 3) * n_0 * n_0 * tau / n_coll_nodes / pow(2, 2.5);

    CollisionNodes collision_nodes(n_coll_nodes, &velocity_grid, 1);
    collision_nodes.randomizeNodes(korobov_coefs);
    collision_nodes.scaleValues();
    collision_nodes.checkOutOfSphere();
    collision_nodes.recalculateRelVelocities();
    int num_file = 0;
    for (double t = 0; t <= t_end; t += tau) {
        std::vector<int> shuffle = generatePermutation(n_coll_nodes);

        for (int pre_i = 0; pre_i < n_coll_nodes; pre_i++) {
            int coll_ix = shuffle[pre_i];
            if (!collision_nodes.isActive(coll_ix)) 
                continue;
            active_collisions++;
            std::tuple<VectorVelocity, VectorVelocity, VectorVelocity> result = collision_nodes.getVelocitiesForInterp(coll_ix);
            InterpNodes got_nodes = velocity_grid.getInterpNodes(result);
            if (got_nodes.r < 0) 
                continue;
            valid_collisions++;
        
            // double energ = velocity_grid.getEnergy();
            double omega = velocity_grid.calculateOmega(got_nodes);
            velocity_grid.updateDistr(got_nodes, omega, integr_constant);
            // double energ_1 = velocity_grid.getEnergy();
            // if (fabs(energ_1 - energ) > 1e-10){
            //     got_nodes.print();
            //     std::cout << "E " << energ_1 - energ << ' '<< omega << std::endl; 
            // }
        }
        // std::cout << n_coll_nodes << ' ' << active_collisions << ' ' << valid_collisions << std::endl;
        // break;
        velocity_grid.saveToFile("checkpoints_new/log" + std::to_string(num_file) + std::to_string(distr_type) + ".txt");
        velocity_grid.printIntegrals();
        std::cout << t << std::endl;
        num_file++;
    }
    
    std::cout << "test 0" << std::endl;
    return 0;
}