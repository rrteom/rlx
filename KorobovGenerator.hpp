#include <vector>

struct KorobovCoefs {
    unsigned int p, b;
    std::vector<unsigned int> coefs;
};

bool is_prime(unsigned int x);
unsigned int next_prime(unsigned int x);
double H(unsigned int z, unsigned int p, unsigned int dim);
std::vector<unsigned int> getHArgmins(unsigned int p, unsigned int dim);
unsigned int getHArgmin(unsigned int p, unsigned int dim);
std::vector<unsigned int> getCoefs(unsigned int b, unsigned int p, unsigned int dim);
std::vector<unsigned int> generateKorobovCoefs(unsigned int N, unsigned int dim);
void print_vector(std::vector<unsigned int> vec);