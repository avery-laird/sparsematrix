//
// Created by avery on 2021-01-16.
//

#include "optimized_triang_test0.cpp" //optsolve0
//#include "optimized_triang_test1.cpp" //optsolve1
//#include "optimized_triang_test2.cpp" //optsolve2
#include "utils.h"
#include "matrix.h"
#include <chrono>
#include <string>
#include <utility>
#include <iostream>


template <class T>
struct TestRecord {
    std::string a_path;
    std::string b_path;
    std::string name;
    int (*solver)(int, int*, int*, T*, T*&);
};

template <class T>
class Test {
public:
    explicit Test(struct TestRecord<T> record)
    : record(std::move(record)) {}

    std::chrono::duration<double> run() {
        // create A and b
        SparseMatrix<T> A(record.a_path);
        Vector<T> b(record.b_path);
        auto begin = std::chrono::high_resolution_clock::now();
        record.solver(A.n, A.col, A.row, A.m, b.m);
        auto end = std::chrono::high_resolution_clock::now();
        return end-begin;
    }

    struct TestRecord<T> record;
};

int main() {

    struct TestRecord<double> torso = {
            "torso1/torso1.mtx",
            "b_for_torso1.mtx",
            "torso test",
            lsolve
    };

    Test<double> naive(torso);
    auto results = naive.run();
    std::cout << results.count();

}
