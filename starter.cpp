#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <cmath>
//#include "testopt.cpp"
#include "testopt_mini.cpp"

//https://en.wikipedia.org/wiki/Machine_epsilon
double epsilon = 1.11e-16;

/*
* Lower triangular solver Lx=b
* L is stored in the compressed column storage format
* Inputs are:
* n : the matrix dimension
* Lp : the column pointer of L
* Li : the row index of L
* Lx : the values of L
* In/Out:
* x : the right hand-side b at start and the solution x at the end.
*/
int lsolve(int n, int *Lp, int *Li, double *Lx, double *&x) {
    int p, j;
    if (!Lp || !Li || !x) return (0); /* check inputs */
    for (j = 0; j < n; j++) {
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
    return (1);
}


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

// FOR TESTING: sympiler's mm reader
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

private:

};

template<class ValueType>
class Vector {
public:
    explicit Vector(std::string path) {
        bool res = inputRHS(path, n, NNZ, m);
        if (!res) throw std::runtime_error("failed to read vector");
    }

    explicit Vector(int n) : n(n) {
        m = new ValueType[n]();
    }

    int n;
    int NNZ{};
    ValueType *m;
};

template<class ValueType>
bool operator==(Vector<ValueType> const &a, Vector<ValueType> const &b) {
    if (a.n != b.n) return false;
    for (int i = 0; i < a.n; i++) {
        if (a.m[i] == b.m[i] ||
            std::isnan(a.m[i]) && std::isnan(b.m[i]))
            continue;
        else return false;
    }
    return true;
}

bool saveVecToMtx(std::string filePath, int n, int nnz, double *&x) {
    std::ofstream output(filePath);
    output << "%%MatrixMarket matrix coordinate real general" << std::endl;
    // write shape
    output << n << " " << 1 << " " << nnz << std::endl;
    for (int i = 0, j = 0; i < n; i++)
        if (x[i] != 0 && !std::isnan(x[i]) && x[i] > epsilon && !std::isinf(x[i])) {
            output << i + 1 << " " << 1 << " " << x[i] << std::endl;
        }
    return true;
}


// FOR TESTING: sympiler function to prepare rhs such that A[1] = b
void rhsInit(int n, int *Ap, int *Ai, double *Ax, double *b) {
    /*generating a rhs that produces a result of all 1 vector*/
    for (int j = 0; j < n; ++j) {
        b[j] = 0;
    }
    for (int c = 0; c < n; ++c) {
        for (int cc = Ap[c]; cc < Ap[c + 1]; ++cc) {
            b[Ai[cc]] += Ax[cc];
        }
    }
}

int main(int argc, char *argv[]) {
    //int n, *Lp = {}, *Li = {};
    //double *Lx = {}, *x = {};
    //inputMatrix("/home/avery/Projects/sparse_matrix_opt/minitest.mtx", n, Lp, Li, Lx);
    SparseMatrix<double> A("/home/avery/Projects/sparse_matrix_opt/rset_example.mtx");
    //SparseMatrix<double> A("/home/avery/Projects/sparse_matrix_opt/torso1/torso1.mtx");

    std::cout << "loaded A..." << std::endl;

    Vector<double> x("/home/avery/Projects/sparse_matrix_opt/rset_example_b.mtx");
    Vector<double> y("/home/avery/Projects/sparse_matrix_opt/rset_example_b.mtx");
    //Vector<double> x("/home/avery/Projects/sparse_matrix_opt/b_for_torso1.mtx");
    //Vector<double> y("/home/avery/Projects/sparse_matrix_opt/b_for_torso1.mtx");

    //Vector<double> x(A.n);
    //rhsInit(A.n, A.col, A.row, A.m, x.m);
    //inputRHS("/home/avery/Projects/sparse_matrix_opt/b_for_torso1.mtx", nb, numNonZero, x);
    std::cout << "loaded b..." << std::endl;
    lsolve(A.n, A.col, A.row, A.m, x.m);
    optsolve(A.n, A.col, A.row, A.m, y.m);
    if (!(y == x)) std::cout << "don't match!" << std::endl;
    else std::cout << "match" << std::endl;
    saveVecToMtx("/home/avery/Projects/sparse_matrix_opt/vol_output_x.mtx", x.n, x.NNZ, x.m);
    saveVecToMtx("/home/avery/Projects/sparse_matrix_opt/vol_output_y.mtx", y.n, y.NNZ, y.m);

    // try in reverse
    //Vector<double> b(x.n);
    // assert that A*x = b
    //spmv_csc(A.n, A.col, A.row, A.m, x.m, b.m);

}