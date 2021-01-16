//
// Created by avery on 2021-01-14.
//

#include <iostream>
#include <Eigen/Core>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/SparseLU>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>


int main(int argc, char *argv[]) {
    typedef Eigen::SparseMatrix<double, Eigen::Lower>SMatrix;
    SMatrix A;
    Eigen::loadMarket(A, "/home/avery/Projects/sparse_matrix_opt/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx");
    std::cout << "load A" << std::endl;
    Eigen::VectorXd b, x, x2;
    Eigen::loadMarketVector(b, "/home/avery/Projects/sparse_matrix_opt/b_for_TSOPF_RS_b678_c2_b.mtx");
    std::cout << "load B" << std::endl;
    Eigen::LeastSquaresConjugateGradient<SMatrix> solver;
    solver.compute(A);
    std::cout << "factorized" << std::endl;
    //if (solver.info() != Eigen::Success)
    //    return 0;
    x = solver.solve(b);
    std::cout << "solved" << std::endl;
    //if (solver.info() != Eigen::Success)
    //    return 0;
    //Eigen::loadMarketVector(x2, "/home/avery/Projects/sparse_matrix_opt/testx.mtx");
    Eigen::saveMarketVector(x2, "/home/avery/Projects/sparse_matrix_opt/eigen_testx.mtx");
    std::cout << "success" << std::endl;
}