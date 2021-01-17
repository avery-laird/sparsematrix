//
// Created by avery on 2021-01-16.
//

#include <iostream>
#include <utility>

class OptimizedSolver {
public:
    OptimizedSolver(std::string A, std::string b)
            : A(std::move(A)), b(std::move(b)) {}

    void generate() {

    }

    std::string A, b;
};


int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "wrong number of arguments" << std::endl;
        return 1;
    }

    /*
     * Get input A matrix and b vector from command line:
     */
    std::string a_path(argv[1]);
    std::string b_path(argv[2]);


}