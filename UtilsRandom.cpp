#include "UtilsRandom.hpp"

std::vector<int> generatePermutation(int permutation_size) {
    std::vector<int> initial, result;
    for (int i = 0; i < permutation_size; i++)
        initial.push_back(i);
    for (int i = 0; i < permutation_size; i++) {
        int rand_index = rand() % initial.size();
        result.push_back(initial[rand_index]);
        std::swap(initial[rand_index], initial[initial.size() - 1]);
        initial.pop_back();
    }
    return result;
}