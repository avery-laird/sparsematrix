#include <fstream>
#include <sstream>
#include <iostream>
#include <set>

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

bool inputRHS(std::string rhsPath, int &n, double *&x) {
    std::ifstream matrixFile(rhsPath);
    std::string line;
    int numberRows, numberCols, numberNonZero;
    while (std::getline(matrixFile, line)) {
        std::istringstream inputstream(line);
        if (line[0] == '%') continue;
        if (!(inputstream >> numberRows >> numberCols >> numberNonZero))
            return false;
        else break;
    }
    x = new double[numberRows]();
    int row, col;
    double val;
    while (std::getline(matrixFile, line)) {
        std::istringstream inputstream(line);
        if (!(inputstream >> row >> col >> val)) return false;
        x[row-1] = val;
    }
    return true;
}

int main(int argc, char *argv[]) {
    int n, *Lp = {}, *Li = {};
    double *Lx = {}, *x = {};
    //inputMatrix("/home/avery/Projects/sparse_matrix_opt/minitest.mtx", n, Lp, Li, Lx);
    inputMatrix("/home/avery/Projects/sparse_matrix_opt/torso1/torso1.mtx", n, Lp, Li, Lx);
    std::cout << "loaded A..." << std::endl;
    int nb, *xCols = {}, *xRows = {};
    inputRHS("/home/avery/Projects/sparse_matrix_opt/b_for_torso1.mtx", nb, x);
    std::cout << "loaded b..." << std::endl;
    lsolve(n, Lp, Li, Lx, x);
    for (int i=0; i < n; i++)
        if (x[i] != 0) {
            std::cout << x[i] << std::endl;
        }

    delete[]Lp;
    delete[]Li;
    delete[]Lx;
    delete []x;
    delete []xCols;
    delete []xRows;
}