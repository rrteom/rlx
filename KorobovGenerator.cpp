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
            prod *= pow(1 - 2 * ((double)mod_buf / p), 2);
        }
        sum += prod;
    }
    return sum;
}

std::vector<unsigned int> getHArgmins(unsigned int p, unsigned int dim) {
    std::vector<unsigned int> argmins = {1};
    double min = H(1, p, dim);
    for (unsigned int z = 2; z <= (p - 1) / 2; z++) {
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

void print_vector(std::vector<unsigned int> vec) {
    for (auto i : vec) 
        std::cout << i << ' ';
    std::cout << std::endl;
    return;
}

int main() {
    unsigned int p_in = 50000;
    std::cout << "p_in init" << std::endl;
    unsigned int p = next_prime(p_in);
    std::cout << "p = " << p << std::endl;
    std::vector<unsigned int> argmins = getHArgmins(p, 8);
    std::cout << "got argmins" << argmins.size() << std::endl;
    for (auto b : argmins) {
        std::cout << "b = " << b << std::endl;
        print_vector(getCoefs(b, p, 8));
    }
    return 0;
}