#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <cmath>
#include <limits>
#include <chrono>
#include "testopt.cpp"
#include "utils.h"
#include "matrix.h"
//#include "testopt_mini.cpp"




/*
* Sparse matrix-vector multiply: y = A*x
* A is stored in the compressed column storage format
* Inputs:
* Ap : the column pointer of A
* Ai : the row index of A
* Ax : the values of A
* x : is a dense vector
* Output:
* y : is a dense vector that stores the result of multiplication
*/
int spmv_csc(int n, int *Ap, int *Ai, double *Ax,
             double *x, double *y) {
    int p, j;
    if (!Ap || !x || !y) return (0);
/* check inputs */
    for (j = 0; j < n; j++) {
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            y[Ai[p]] += Ax[p] * x[j];
        }
    }
    return (1);
}


bool saveVecToMtx(std::string filePath, int n, int nnz, double *&x) {
    std::ofstream output(filePath);
    output << "%%MatrixMarket matrix coordinate real general" << std::endl;
    // write shape
    output << n << " " << 1 << " " << nnz << std::endl;
    for (int i = 0, j = 0; i < n; i++)
        if (x[i] != 0 && !std::isnan(x[i]) && x[i] > std::numeric_limits<double>::min() && !std::isinf(x[i])) {
            output << i + 1 << " " << 1 << " " << x[i] << std::endl;
        }
    return true;
}


int main(int argc, char *argv[]) {
    //int n, *Lp = {}, *Li = {};
    //double *Lx = {}, *x = {};
    //inputMatrix("/home/avery/Projects/sparse_matrix_opt/minitest.mtx", n, Lp, Li, Lx);
    //SparseMatrix<double> A("/home/avery/Projects/sparse_matrix_opt/rset_example.mtx");
    SparseMatrix<double> A("/home/avery/Projects/sparse_matrix_opt/torso1/torso1.mtx");

    std::cout << "loaded A..." << std::endl;

    //Vector<double> x("/home/avery/Projects/sparse_matrix_opt/rset_example_b.mtx");
    //Vector<double> y("/home/avery/Projects/sparse_matrix_opt/rset_example_b.mtx");
    Vector<double> x("/home/avery/Projects/sparse_matrix_opt/b_for_torso1.mtx");
    Vector<double> y("/home/avery/Projects/sparse_matrix_opt/b_for_torso1.mtx");

    //Vector<double> x(A.n);
    //rhsInit(A.n, A.col, A.row, A.m, x.m);
    //inputRHS("/home/avery/Projects/sparse_matrix_opt/b_for_torso1.mtx", nb, numNonZero, x);
    std::cout << "loaded b..." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();
    lsolve(A.n, A.col, A.row, A.m, x.m);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_lsolve = end-begin;
    std::cout << "lsolve elapsed time: " << elapsed_lsolve.count() << std::endl;

    begin = std::chrono::high_resolution_clock::now();
    optsolve(A.n, A.col, A.row, A.m, y.m);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_optsolve = end-begin;
    std::cout << "optsolve elapsed time: " << elapsed_optsolve.count() << std::endl;

    if (!(y == x)) std::cout << "don't match!" << std::endl;
    else std::cout << "match" << std::endl;
    saveVecToMtx("/home/avery/Projects/sparse_matrix_opt/vol_output_x.mtx", x.n, x.NNZ, x.m);
    saveVecToMtx("/home/avery/Projects/sparse_matrix_opt/vol_output_y.mtx", y.n, y.NNZ, y.m);

    // try in reverse
    //Vector<double> b(x.n);
    // assert that A*x = b
    //spmv_csc(A.n, A.col, A.row, A.m, x.m, b.m);

}