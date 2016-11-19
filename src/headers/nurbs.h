
#ifndef NURBS_H
#define NURBS_H

#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>

namespace nurbs{

typedef Eigen::Triplet<double> trip;
typedef Eigen::SparseMatrix<double> spMat;

struct NurbsBase2D
{
    //
    NurbsBase2D(Eigen::VectorXi u_knots, Eigen::VectorXi v_knots,
              Eigen::VectorXd weights,
              int degree_u=3, int degree_v=3);
    int degree_u = degree_u;
    int degree_v = degree_v;
    Eigen::VectorXi u_knots;
    Eigen::VectorXi v_knots;
    Eigen::VectorXd weights;

    std::vector<std::function<double(double)>> u_functions;
    std::vector<std::function<double(double)>> v_functions;

    std::vector<std::function<double(double)>> Du_functions;
    std::vector<std::function<double(double)>> Dv_functions;

    std::vector<std::function<double(double)>> DDu_functions;
    std::vector<std::function<double(double)>> DDv_functions;

    void computeFirstDerivatives();
    void computeSecondDerivatives();

    Eigen::VectorXd getInfluenceVector(Eigen::Vector2d u);
    spMat getInfluenceMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U);

    Eigen::VectorXd getDuVector(Eigen::Vector2d u);
    spMat getDuMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U);

    Eigen::VectorXd getDvVector(Eigen::Vector2d u);
    spMat getDvMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U);
};

struct NurbsBase1D
{
    NurbsBase1D(Eigen::VectorXi u_knots, Eigen::VectorXd weights, int degree_u=3);
    int degree_u = degree_u;
    Eigen::VectorXi u_knots;
    Eigen::VectorXd weights;
    std::vector<std::function<double(double)>> u_functions;
    std::vector<std::function<double(double)>> Du_functions;
    std::vector<std::function<double(double)>> DDu_functions;

    void computeFirstDerivatives();
    void computeSecondDerivatives();

    Eigen::VectorXd getInfluenceVector(double u);
    spMat getInfluenceMatrix(Eigen::VectorXd u);

    Eigen::VectorXd getDuVector(double u);
    spMat getDuMatrix(Eigen::VectorXd u);

};

}

#endif