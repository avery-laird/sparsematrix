//
// Created by avery on 2021-01-16.
//

#ifndef SPARSE_MATRIX_OPT_MATRIX_H
#define SPARSE_MATRIX_OPT_MATRIX_H

#include <string>
#include <stdexcept>
#include "utils.h"

template<class ValueType>
class Vector {
public:
    explicit Vector(std::string path) {
        bool res = inputRHS(path, n, NNZ, m);
        if (!res) throw std::runtime_error("failed to read vector");
    }
    explicit Vector(int n) {
        m = new ValueType[n]();
    }
    int n;
    int NNZ;
    ValueType *m;
};


template<class ValueType>
class SparseMatrix {
public:
    explicit SparseMatrix(std::string path) {
        bool res = readMatrix(path, n, NNZ, col, row, m);
        if (!res) throw std::runtime_error("failed to read matrix");
    }

    int *row{}; // row index
    int *col{}; // column pointers
    int n;
    int NNZ;
    ValueType *m;
};


template<class ValueType>
bool operator==(Vector<ValueType> const &a, Vector<ValueType> const &b);

#endif //SPARSE_MATRIX_OPT_MATRIX_H
