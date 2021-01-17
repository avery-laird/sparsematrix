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

struct element {
    int row;
    int col;
    double val;

    element(int row, int col, double val)
            : row(row), col(col), val(val) {}
};

/*
 * Given a filepath, copy the non-zero entries of the matrix into memory.
 * The matrix must be in matrix market format.
 * The matrix is stored in memory in CSC format.
 * The file preamble contains lines starting with %, which are ignored
 * The first line which does not start with % is the form: <# rows> <# cols> <nnz>
 * All subsequent lines are of the form <row> <col> <val>
 *
 */
bool inputMatrix(std::string filePath, int &n, int *&Lp, int *&Li, double *&Lx) {
    // TODO: ignore everything above diagonal
    std::ifstream matrixFile(filePath);
    std::string line;
    int numberRows, numberCols, numberNonZero;
    while (std::getline(matrixFile, line)) {
        std::istringstream inputstream(line);
        if (line[0] == '%') continue;
        if (!(inputstream >> numberRows >> numberCols >> numberNonZero))
            return false;
        else break;
    }
    // assume L is square
    //if (numberRows != numberCols) return false; // TODO: re-enable
    // allocate memory for L and b
    Lp = new int[numberCols]; // column pointers
    Li = new int[numberNonZero]; // row pointers
    Lx = new double[numberNonZero]; // values
    n = numberCols;
    // 2nd phase: read the actual matrix.
    // Imagine the full matrix is stretched out and
    // arranged left-to-right by column. First store
    // the start of each non-zero stretch in Lp.
    // Then, remove all the zeros, and store nz elements
    // in Lx. Finally, for each element in Lx, store the
    // row in Li.
    // The matrixmarket format doesn't guarantee that elements
    // are in any strictly increasing order.
    int indexp = 0;
    int row, col;
    double val;
    // first, read all values.
    // https://stackoverflow.com/questions/2620862/using-custom-stdset-comparator
    auto cmp = [numberCols](element a, element b) { return a.col * numberCols + a.row < b.col * numberCols + b.row; };
    std::multiset<element, decltype(cmp)> entries(cmp);
    while (std::getline(matrixFile, line)) {
        std::istringstream inputstream(line);
        if (!(inputstream >> row >> col >> val)) return false;
        if (col > row) continue;
        entries.emplace(row - 1, col - 1, val);
    }
    int currentCol = 0;
    int countNonZero = 0;
    int i = 0;
    Lp[0] = 0;
    int lpp = 0;
    for (auto entry : entries) {
        Lx[i] = entry.val; // always store the value no matter what
        // next, store the row
        Li[i++] = entry.row;
        // Lp is groups of two ranges, start col and end col
        // every time the column increases, add the index of
        // the last entry in Lx
        for (; currentCol < entry.col; currentCol++, countNonZero = 0) {
            Lp[lpp + 1] = Lp[lpp] + countNonZero;
            lpp++;
        }
        countNonZero++;
    }
    return true;
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