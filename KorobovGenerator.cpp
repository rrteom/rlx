#include <iostream>
#include <vector>
#include <cmath>

bool is_prime(unsigned int x) {
    for (unsigned int i = 2; i <= x / 2; i++)
        if (x % i == 0)
            return false;
    return true;
}

unsigned int next_prime(unsigned int x) {
    unsigned int ans = x;
    while (!is_prime(ans)) ans++;
    return ans;
}

double H(unsigned int z, unsigned int p, unsigned int dim) {
    double sum = 0;
    for (unsigned int k = 1; k <= (p - 1) / 2; k++) {
        unsigned int mod_buf = k;
        double prod = 1;
        for (unsigned int n = 1; n <= dim; n++) {
            mod_buf = mod_buf * z % p;
            prod *= 1 - 2 * ((double)mod_buf / p);
        }
        sum += prod * prod;
    }
    return sum;
}

std::vector<unsigned int> getHArgmins(unsigned int p, unsigned int dim) {
    std::vector<unsigned int> argmins = {1};
    double min = H(1, p, dim);
    for (unsigned int z = 2; z <= (p - 1) / 2; z++) {
        if (z % 1000 == 0)
            std::cout << "z = " << z << std::endl;
        if (H(z, p, dim) == min)
            argmins.push_back(z);
        else if (H(z, p, dim) < min) {
            min = H(z, p, dim);
            argmins.clear();
            argmins.push_back(z);
        }
    }
    return argmins;
}

unsigned int getHArgmin(unsigned int p, unsigned int dim) {
    unsigned int argmin = 1;
    double min = H(1, p, dim);
    for (unsigned int z = 2; z <= (p - 1) / 2; z++) {
        if (z % 1000 == 0)
            std::cout << "z = " << z << std::endl;
        if (H(z, p, dim) < min) {
            min = H(z, p, dim);
            argmin = z;
        }
    }
    return argmin;
}

std::vector<unsigned int> getCoefs(unsigned int b, unsigned int p, unsigned int dim) {
    // coefs for Korobov generator
    std::vector<unsigned int> coefs = {1, b};
    unsigned int mod_buf = b;
    for (unsigned int i = 2; i < dim; i++) {
        mod_buf = b * mod_buf % p;
        coefs.push_back(mod_buf);
    }
    return coefs;
}

std::vector<unsigned int> generateKorobovCoefs(unsigned int p_in, unsigned int dim, unsigned int n_nodes) {
    unsigned int p = next_prime(p_in);
    unsigned int b = getHArgmin(p, dim);
    return getCoefs(b, p, dim);
}

void print_vector(std::vector<unsigned int> vec) {
    for (auto i : vec) 
        std::cout << i << ' ';
    std::cout << std::endl;
    return;
}