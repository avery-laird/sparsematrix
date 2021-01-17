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
            : record(record), A(record.a_path), b(record.b_path), x(A.n) {}

    std::chrono::duration<double> run(int (*solver)(int, int *, int *, T *, T *&)) {
        auto begin = std::chrono::high_resolution_clock::now();
        solver(A.n, A.col, A.row, A.m, x.m);
        auto end = std::chrono::high_resolution_clock::now();
        return end - begin;
    }

    struct TestRecord<T> record;
    SparseMatrix<T> A;
    Vector<T> b;
    Vector<T> x;
};


int main() {
    /*
     * Perform benchmarks, and verify that solutions
     * agree with the naive implementation
     */
    Test<double> torso({
                               "../torso1/torso1.mtx",
                               "../b_for_torso1.mtx",
                               "torso test"
                       });
    auto results = torso.run(lsolve);
    Vector<double> solution1 = torso.x;
    std::cout << "lsolve, torso1: " << results.count() << std::endl;

    results = torso.run(optsolve0);
    Vector<double> solution2 = torso.x;
    std::cout << "optsolve, torso1: " << results.count() << std::endl;

    if (solution1 == solution2)
        std::cout << "match!";
    else std::cout << "don't match";
    std::cout << std::endl;

    /**
     * Begin test of the second matrix
     */
    Test<double> tsopf({
        "../TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx",
        "../b_for_TSOPF_RS_b678_c2_b.mtx"
    });
    results = tsopf.run(lsolve);
    solution1 = tsopf.x;
    std::cout << "lsolve, TSOPF: " << results.count() << std::endl;

    results = tsopf.run(optsolve1);
    solution2 = tsopf.x;
    std::cout << "optsolve1, TSOPF: " << results.count() << std::endl;

    if (solution1 == solution2)
        std::cout << "match!";
    else std::cout << "don't match";
    std::cout << std::endl;


}
