//
// Created by avery on 2021-01-16.
//

#include "optimized_triang_test0.cpp" //optsolve0=torso
#include "optimized_triang_test1.cpp" //optsolve1=TSOPF_RS_b678_c2
#include "optimized_triang_test2.cpp" //optsolve2=example
#include "utils.h"
#include "matrix.h"
#include <chrono>
#include <string>
#include <utility>
#include <iostream>


template<class T>
struct TestRecord {
    std::string a_path;
    std::string b_path;
    std::string name;
};

template<class T>
class Test {
public:
    explicit Test(struct TestRecord<T> record)
            : record(record), A(record.a_path), b(record.b_path) {}

    std::pair<Vector<T>, std::chrono::duration<double>> run(int (*solver)(int, int *, int *, T *, T *&)) {
        // allocate a new vector for the solution
        Vector<T> x(A.n);
        auto begin = std::chrono::high_resolution_clock::now();
        solver(A.n, A.col, A.row, A.m, x.m);
        auto end = std::chrono::high_resolution_clock::now();
        return {x, end - begin};
    }

    struct TestRecord<T> record;
    SparseMatrix<T> A;
    Vector<T> b;
};

int main() {


    Test<double> torso({
                               "../torso1/torso1.mtx",
                               "../b_for_torso1.mtx",
                               "torso test"
                       });
    auto results = torso.run(lsolve);
    Vector<double> solution1 = results.first;
    std::cout << results.second.count() << std::endl;

    results = torso.run(optsolve0);
    Vector<double> solution2 = results.first;
    std::cout << results.second.count() << std::endl;

    if (solution1 == solution2)
        std::cout << "match!";
    else std::cout << "don't match";
    std::cout << std::endl;

}
