//
// Created by avery on 2021-01-16.
//

#ifndef SPARSE_MATRIX_OPT_UTILS_H
#define SPARSE_MATRIX_OPT_UTILS_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

int lsolve(int n, int *Lp, int *Li, double *Lx, double *&x);


// FOR TESTING: sympiler's mm reader
// TODO replace with my own later
/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
template<class ValueType>
bool readMatrix(std::string fName, int &n, int &NNZ, int *&col,
                int *&row, ValueType *&val) {
    /*This function reads the input matrix from "fName" file and
     * allocate memory for matrix A, L and U.
     * - The input file is a coordinate version and e
     * ach row of the file shows (col, row, nnz)
     * - The matrices are zero-indexed
     */

    std::ifstream inFile;
    inFile.open(fName);
    std::string line, banner, mtx, crd, arith, sym;
    /*  File format:
     *    %%MatrixMarket matrix coordinate real general/symmetric/...
     *    % ...
     *    % (optional comments)
     *    % ...
     *    #rows    #non-zero
     *    Triplet in the rest of lines: row    col    value
     */
    std::getline(inFile, line);
    for (unsigned i = 0; i < line.length(); line[i] = tolower(line[i]), i++);
    std::istringstream iss(line);
    if (!(iss >> banner >> mtx >> crd >> arith >> sym)) {
        std::cout << "Invalid header (first line does not contain 5 tokens)\n";
        return false;
    }

    if (banner.compare("%%matrixmarket")) {
        std::cout << "Invalid header (first token is not \"%%%%MatrixMarket\")\n";
        return false;
    }
    if (mtx.compare("matrix")) {
        std::cout << "Not a matrix; this driver cannot handle that.\"\n";
        return false;
    }
    if (crd.compare("coordinate")) {
        std::cout << "Not in coordinate format; this driver cannot handle that.\"\n";
        return false;
    }
    if (arith.compare("real")) {
        if (!arith.compare("complex")) {
            std::cout << "Complex matrix; use zreadMM instead!\n";
            return false;
        } else if (!arith.compare("pattern")) {
            std::cout << "Pattern matrix; values are needed!\n";
            return false;
        } else {
            std::cout << "Unknown arithmetic\n";
            return false;
        }
    }
    while (!line.compare(0, 1, "%")) {
        std::getline(inFile, line);
    }
    std::istringstream issDim(line);
    if (!(issDim >> n >> n >> NNZ)) {
        std::cout << "The matrix dimension is missing\n";
        return false;
    }
    if (n <= 0 || NNZ <= 0)
        return false;
    col = new int[n + 1]();
    // colL = new int[n + 1]; colU = new int[n + 1];
    row = new int[NNZ];
    // rowL = new int[factorSize]; rowU = new int[factorSize];
    val = new ValueType[NNZ];
    // valL = new double[factorSize]; valU = new double[factorSize];
    if (!val || !col || !row)
        return false;
    //Initializing the result vector
    int y, x, colCnt = 0, nnzCnt = 0;
    double value;

    col[0] = 0;
    for (int i = 0; nnzCnt < NNZ;) {//Reading from file row by row
        inFile >> x;
        x--;
        inFile >> y;
        y--;//zero indexing
        inFile >> value;
        if (y > n)
            return false;
        if (y == i) {
            val[nnzCnt] = value;
            row[nnzCnt] = x;
            colCnt++;
            nnzCnt++;
        } else {//New col
            col[i + 1] = col[i] + colCnt;
            i++;//next iteration
            colCnt = 1;
            val[nnzCnt] = value;
            row[nnzCnt] = x;
            nnzCnt++;
        }

    }
    col[n] = col[n - 1] + colCnt;//last col

    return true;
}

bool nearlyEqual(float a, float b);

template<class ValueType>
bool inputRHS(std::string rhsPath, int &n, int &numberNonZero, ValueType *&x) {
    std::ifstream matrixFile(rhsPath);
    std::string line;
    int numberCols;
    while (std::getline(matrixFile, line)) {
        std::istringstream inputstream(line);
        if (line[0] == '%') continue;
        if (!(inputstream >> n >> numberCols >> numberNonZero))
            return false;
        else break;
    }
    x = new ValueType[n]();
    int row, col;
    double val;
    while (std::getline(matrixFile, line)) {
        std::istringstream inputstream(line);
        if (!(inputstream >> row >> col >> val)) return false;
        x[row - 1] = val;
    }
    return true;
}

#endif //SPARSE_MATRIX_OPT_UTILS_H
